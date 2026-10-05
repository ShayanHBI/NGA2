"""
Reference two-phase relaxation: p_relax, pT_relax, pTg_relax.

Each entry point is a self-contained nonlinear solve on the cell's extensive
state.  No linearization, no basis decomposition, no step ordering: a relax
routine is defined by the equilibrium conditions it imposes and by the cell
quantities it conserves, and it solves for every coordinate those conditions
leave free.

    p_relax    p_l - p_g = p_jump                        free: V_l  (E_l slaved,
                                                                no heat allowed)
    pT_relax   p_l - p_g = p_jump,  T_l = T_g            free: V_l, E_l
    pTg_relax  p_l - p_g = p_jump,  T_l = T_g,  g_l = mu_v
                                                         free: V_l, E_l, m_l

Conserved in all three: V_l + V_g, E_l + E_g + p_jump*V_l, m_l + m_v, m_a.

Liquid: NASG.  Gas: ideal mixture of vapour and air.
Parameters default to the Sembian blastwave NASG case.
"""

import numpy as np
from scipy.optimize import root, brentq
from scipy.integrate import solve_ivp


# ----------------------------------------------------------------------
# equations of state
# ----------------------------------------------------------------------

class NASG:
    """Noble-Abel stiffened gas, single component."""

    def __init__(self, gamma=1.19, cv=3610.0, pinf=7.028e8, b=6.61e-4,
                 q=-1177788.0, qp=0.0):
        self.gamma, self.cv, self.pinf, self.b, self.q, self.qp = \
            gamma, cv, pinf, b, q, qp
        self.cp = gamma * cv
        self.R = (gamma - 1.0) * cv

    def T(self, V, E, m):
        v, e = V / m, E / m
        return (e - self.q - self.pinf * (v - self.b)) / self.cv

    def p(self, V, E, m):
        v = V / m
        return (self.gamma - 1.0) * self.cv * self.T(V, E, m) / (v - self.b) \
            - self.pinf

    def S(self, V, E, m):
        """Extensive entropy.

        The reference constant MUST match g_Tp below.  Writing this as
        cv*ln(e-q-pinf*(v-b)) + R*ln(v-b) + qp -- the form usually quoted --
        differs from the (T,p) form by cv*ln(cv) + R*ln(R), which for the
        Sembian liquid is 3.4e4 J/(kg K), i.e. a 1.5e7 J/kg offset in g_l at
        430 K.  That dwarfs the latent heat and destroys chemical equilibrium.
        """
        v, e = V / m, E / m
        cvT = e - self.q - self.pinf * (v - self.b)
        vb = v - self.b
        if cvT <= 0.0 or vb <= 0.0:
            return -np.inf
        return m * (self.cv * np.log(cvT / self.cv) + self.R * np.log(vb)
                    - self.R * np.log(self.R) + self.qp)

    def g_Tp(self, T, p):
        """Specific Gibbs energy from (T, p).  Matches nasg_get_g_from_p_T."""
        h = self.cp * T + self.b * p + self.q
        s = self.cp * np.log(T) - (self.gamma - 1.0) * self.cv \
            * np.log(p + self.pinf) + self.qp
        return h - T * s

    def rho_Tp(self, T, p):
        return 1.0 / (self.b + (self.gamma - 1.0) * self.cv * T / (p + self.pinf))

    def E_of(self, T, m, V):
        v = V / m
        return m * (self.q + self.cv * T + self.pinf * (v - self.b))


class IdealSpecies:
    def __init__(self, gamma, cv, q=0.0, qp=0.0):
        self.gamma, self.cv, self.q, self.qp = gamma, cv, q, qp
        self.cp = gamma * cv
        self.R = (gamma - 1.0) * cv

    def h(self, T):
        return self.cp * T + self.q

    def s(self, T, p_partial):
        return self.cp * np.log(T) - self.R * np.log(p_partial) + self.qp

    def g(self, T, p_partial):
        """Pure-species Gibbs energy at whatever pressure it is evaluated at."""
        return self.h(T) - T * self.s(T, p_partial)


class IdealMixture:
    """Ideal mixture of ideal gases.  Species 0 is the condensable (vapour)."""

    def __init__(self, species):
        self.sp = species

    def T(self, E, masses):
        num = E - sum(m * s.q for m, s in zip(masses, self.sp))
        den = sum(m * s.cv for m, s in zip(masses, self.sp))
        return num / den

    def N(self, masses):
        return sum(m * s.R for m, s in zip(masses, self.sp))

    def p(self, V, E, masses):
        return self.N(masses) * self.T(E, masses) / V

    def x(self, masses):
        N = self.N(masses)
        return np.array([m * s.R / N for m, s in zip(masses, self.sp)])

    def S(self, V, E, masses):
        T = self.T(E, masses)
        if T <= 0.0 or V <= 0.0:
            return -np.inf
        pg = self.p(V, E, masses)
        xs = self.x(masses)
        tot = 0.0
        for m, s, xi in zip(masses, self.sp, xs):
            if m <= 0.0:
                continue
            if xi <= 0.0:
                return -np.inf
            tot += m * s.s(T, xi * pg)
        return tot

    def mu(self, i, V, E, masses):
        """Mass-specific chemical potential of species i in the mixture."""
        T = self.T(E, masses)
        pg = self.p(V, E, masses)
        p_i = self.x(masses)[i] * pg
        return self.sp[i].g(T, p_i)

    def E_of(self, T, masses):
        return sum(m * (s.cv * T + s.q) for m, s in zip(masses, self.sp))


# ----------------------------------------------------------------------
# cell state
# ----------------------------------------------------------------------

class Cell:
    """A mixed cell.  State: V_l, E_l, m_l, m_v, m_a (V, E totals conserved)."""

    def __init__(self, liq, gas, V, E, V_l, E_l, m_l, m_v, m_a, p_jump=0.0):
        self.liq, self.gas = liq, gas
        self.V, self.E = V, E          # totals (E includes interface energy)
        self.V_l, self.E_l = V_l, E_l
        self.m_l, self.m_v, self.m_a = m_l, m_v, m_a
        self.p_jump = p_jump

    @classmethod
    def from_primitive(cls, liq, gas, V, VF, T_l, p_l, T_g, p_g, x_v,
                       p_jump=0.0):
        V_l, V_g = VF * V, (1.0 - VF) * V
        m_l = liq.rho_Tp(T_l, p_l) * V_l
        E_l = liq.E_of(T_l, m_l, V_l)
        Ntot = p_g * V_g / T_g
        m_v = x_v * Ntot / gas.sp[0].R
        m_a = (1.0 - x_v) * Ntot / gas.sp[1].R
        E_g = gas.E_of(T_g, [m_v, m_a])
        E = E_l + E_g + p_jump * V_l
        return cls(liq, gas, V, E, V_l, E_l, m_l, m_v, m_a, p_jump)

    # -- derived -------------------------------------------------------
    def V_g(self, V_l=None):
        return self.V - (self.V_l if V_l is None else V_l)

    def E_g(self, V_l=None, E_l=None):
        V_l = self.V_l if V_l is None else V_l
        E_l = self.E_l if E_l is None else E_l
        return self.E - E_l - self.p_jump * V_l

    def masses_g(self, m_l=None):
        m_l = self.m_l if m_l is None else m_l
        return [self.m_w - m_l, self.m_a]

    @property
    def m_w(self):
        return self.m_l + self.m_v

    def state(self, V_l=None, E_l=None, m_l=None):
        V_l = self.V_l if V_l is None else V_l
        E_l = self.E_l if E_l is None else E_l
        m_l = self.m_l if m_l is None else m_l
        return V_l, E_l, m_l, self.V - V_l, self.E - E_l - self.p_jump * V_l, \
            [self.m_w - m_l, self.m_a]

    def report(self, V_l=None, E_l=None, m_l=None):
        Vl, El, ml, Vg, Eg, mg = self.state(V_l, E_l, m_l)
        Tl = self.liq.T(Vl, El, ml)
        pl = self.liq.p(Vl, El, ml)
        Tg = self.gas.T(Eg, mg)
        pg = self.gas.p(Vg, Eg, mg)
        gl = self.liq.g_Tp(Tl, pl)
        muv = self.gas.mu(0, Vg, Eg, mg)
        return dict(V_l=Vl, E_l=El, m_l=ml, m_v=mg[0], VF=Vl / self.V,
                    T_l=Tl, T_g=Tg, p_l=pl, p_g=pg,
                    x_v=self.gas.x(mg)[0], g_l=gl, mu_v=muv,
                    F_M=pl - pg - self.p_jump, F_T=Tl - Tg, F_C=gl - muv,
                    S=self.liq.S(Vl, El, ml) + self.gas.S(Vg, Eg, mg))

    def copy_with(self, V_l, E_l, m_l):
        m_v = self.m_w - m_l
        return Cell(self.liq, self.gas, self.V, self.E, V_l, E_l, m_l, m_v,
                    self.m_a, self.p_jump)


# ----------------------------------------------------------------------
# the three relaxations
# ----------------------------------------------------------------------

def p_relax(cell, pI_model="gas", rtol=1e-12):
    """Mechanical relaxation.

    Volume crosses the interface; no heat and no mass.  The no-heat constraint
    slaves the energy to the volume through the interface work,

        dE_l = -p_I dV_l,   dE_g = +(p_I - p_jump) dV_l,

    so there is one free coordinate, V_l, and one end condition,
    p_l - p_g = p_jump.  p_I is a closure; the answer is insensitive to it at
    first order in the force (see pI_model comparison in the self-test).
    """
    Vl0, El0, ml = cell.V_l, cell.E_l, cell.m_l
    mg = [cell.m_v, cell.m_a]

    def pI_of(Vl, El):
        pl = cell.liq.p(Vl, El, ml)
        pg = cell.gas.p(cell.V - Vl, cell.E - El - cell.p_jump * Vl, mg)
        if pI_model == "liquid":
            return pl
        if pI_model == "gas":
            return pg + cell.p_jump
        if pI_model == "mean":
            return 0.5 * (pl + pg + cell.p_jump)
        raise ValueError(pI_model)

    def rhs(Vl, y):
        return [-pI_of(Vl, y[0])]

    def El_at(Vl):
        if abs(Vl - Vl0) < 1e-300:
            return El0
        sol = solve_ivp(rhs, (Vl0, Vl), [El0], rtol=1e-12, atol=1e-10,
                        dense_output=True, max_step=abs(Vl - Vl0) / 8.0)
        return float(sol.y[0, -1])

    def resid(Vl):
        El = El_at(Vl)
        pl = cell.liq.p(Vl, El, ml)
        pg = cell.gas.p(cell.V - Vl, cell.E - El - cell.p_jump * Vl, mg)
        return pl - pg - cell.p_jump

    lo, hi = _bracket(resid, Vl0, cell.V)
    Vl = brentq(resid, lo, hi, xtol=cell.V * 1e-15, rtol=8.9e-16)
    return cell.copy_with(Vl, El_at(Vl), ml)


def pT_relax(cell):
    """Mechanical + thermal relaxation, solved as one system.

    Imposes BOTH p_l - p_g = p_jump and T_l = T_g.  Two end conditions, two
    free coordinates (V_l, E_l): the volume exchange and the heat exchange are
    solved together, simultaneously, from the cell's conserved totals.  No
    prior p_relax is needed and none is assumed.
    """
    ml = cell.m_l
    mg = [cell.m_v, cell.m_a]

    def resid(u):
        Vl, El = u
        Vg, Eg = cell.V - Vl, cell.E - El - cell.p_jump * Vl
        pl = cell.liq.p(Vl, El, ml)
        pg = cell.gas.p(Vg, Eg, mg)
        Tl = cell.liq.T(Vl, El, ml)
        Tg = cell.gas.T(Eg, mg)
        return [(pl - pg - cell.p_jump) / max(abs(pl), 1.0),
                (Tl - Tg) / max(abs(Tl), 1.0)]

    sol = _solve(resid, [cell.V_l, cell.E_l], "pT_relax")
    return cell.copy_with(sol[0], sol[1], ml)


def pTg_relax(cell):
    """Full relaxation: mechanical + thermal + chemical, solved as one system.

    Imposes p_l - p_g = p_jump, T_l = T_g and g_l = mu_v.  Three end
    conditions, three free coordinates (V_l, E_l, m_l).
    """
    def resid(u):
        Vl, El, ml = u
        mv = cell.m_w - ml
        Vg, Eg = cell.V - Vl, cell.E - El - cell.p_jump * Vl
        mg = [mv, cell.m_a]
        pl = cell.liq.p(Vl, El, ml)
        pg = cell.gas.p(Vg, Eg, mg)
        Tl = cell.liq.T(Vl, El, ml)
        Tg = cell.gas.T(Eg, mg)
        gl = cell.liq.g_Tp(Tl, pl)
        muv = cell.gas.mu(0, Vg, Eg, mg)
        return [(pl - pg - cell.p_jump) / max(abs(pl), 1.0),
                (Tl - Tg) / max(abs(Tl), 1.0),
                (gl - muv) / max(abs(gl), 1.0)]

    sol = _solve(resid, [cell.V_l, cell.E_l, cell.m_l], "pTg_relax")
    return cell.copy_with(*sol)


def _solve(resid, x0, tag, tries=6):
    """Robust scaled root find: hybr, then Newton polish on the raw residual."""
    x = np.array(x0, float)
    best, best_n = None, np.inf
    for k in range(tries):
        sol = root(resid, x, method="hybr")
        n = np.linalg.norm(resid(sol.x))
        if n < best_n:
            best, best_n = sol.x.copy(), n
        if n < 1e-13:
            return sol.x
        # restart from the best point, nudged, to escape a stalled step
        x = best * (1.0 + 1e-9 * (k + 1))
    # Newton polish with a numerical Jacobian
    x = best.copy()
    for _ in range(60):
        f = np.array(resid(x))
        if np.linalg.norm(f) < 1e-14:
            break
        J = np.zeros((len(x), len(x)))
        for j in range(len(x)):
            h = abs(x[j]) * 1e-8 or 1e-12
            xp = x.copy(); xp[j] += h
            J[:, j] = (np.array(resid(xp)) - f) / h
        try:
            dx = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError:
            break
        lam = 1.0
        for _ in range(40):
            xn = x + lam * dx
            if np.linalg.norm(np.array(resid(xn))) < np.linalg.norm(f):
                x = xn
                break
            lam *= 0.5
        else:
            break
    if np.linalg.norm(np.array(resid(x))) > 1e-9:
        raise RuntimeError(f"{tag} failed, residual "
                           f"{np.linalg.norm(np.array(resid(x))):.3e}")
    return x


def _bracket(f, x0, V):
    """Expand a bracket around x0 inside (0, V) until f changes sign."""
    f0 = f(x0)
    step = V * 1e-6
    for _ in range(80):
        for lo, hi in ((max(x0 - step, V * 1e-12), x0),
                       (x0, min(x0 + step, V * (1 - 1e-12)))):
            try:
                if f(lo) * f(hi) <= 0.0:
                    return lo, hi
            except (ValueError, ZeroDivisionError, FloatingPointError):
                pass
        step *= 2.0
        if step > V:
            break
    raise RuntimeError("could not bracket the mechanical root (f0=%g)" % f0)


# ----------------------------------------------------------------------
# self-test
# ----------------------------------------------------------------------

def _fmt(tag, r):
    return (f"{tag:<12s} VF={r['VF']:.8f}  T_l={r['T_l']:9.4f}  T_g={r['T_g']:9.4f}"
            f"  p_l={r['p_l']:12.4f}  p_g={r['p_g']:12.4f}  Y_v-ish x_v={r['x_v']:.6f}\n"
            f"{'':12s} F_M={r['F_M']: .6e}  F_T={r['F_T']: .6e}  F_C={r['F_C']: .6e}"
            f"   S={r['S']:.10f}")


def main():
    liq = NASG()
    vap = IdealSpecies(1.47, 955.0, 2077616.0, 14317.0)
    air = IdealSpecies(1.40, 718.0, 0.0, 0.0)
    gas = IdealMixture([vap, air])

    cell = Cell.from_primitive(liq, gas, V=1.0e-3, VF=0.4,
                               T_l=430.0, p_l=3.4e5,
                               T_g=400.0, p_g=2.8e5, x_v=0.05)

    r0 = cell.report()
    print("=" * 78)
    print(_fmt("initial", r0))

    c_p = p_relax(cell)
    c_pT = pT_relax(cell)
    c_pTg = pTg_relax(cell)
    print("-" * 78)
    print(_fmt("p_relax", c_p.report()))
    print(_fmt("pT_relax", c_pT.report()))
    print(_fmt("pTg_relax", c_pTg.report()))

    # ---- 1. each routine meets its own end conditions -----------------
    print("=" * 78)
    print("1. end conditions actually met (nonlinear residuals)")
    print(f"   p_relax    |F_M| = {abs(c_p.report()['F_M']):.3e}")
    rT = c_pT.report()
    print(f"   pT_relax   |F_M| = {abs(rT['F_M']):.3e}   |F_T| = {abs(rT['F_T']):.3e}")
    rG = c_pTg.report()
    print(f"   pTg_relax  |F_M| = {abs(rG['F_M']):.3e}   |F_T| = {abs(rG['F_T']):.3e}"
          f"   |F_C| = {abs(rG['F_C']):.3e}")

    # ---- 2. pT_relax is correct on its own ----------------------------
    print("=" * 78)
    print("2. pT_relax does not depend on whether p_relax ran first")
    a = pT_relax(cell).report()
    b = pT_relax(c_p).report()
    print(f"   from raw state      : V_l={a['V_l']:.15e}  E_l={a['E_l']:.15e}")
    print(f"   after p_relax       : V_l={b['V_l']:.15e}  E_l={b['E_l']:.15e}")
    print(f"   relative difference : V_l {abs(a['V_l']/b['V_l']-1):.3e}"
          f"   E_l {abs(a['E_l']/b['E_l']-1):.3e}")

    print("   pTg_relax likewise, from three different starting points")
    g1 = pTg_relax(cell).report()
    g2 = pTg_relax(c_p).report()
    g3 = pTg_relax(c_pT).report()
    for tag, g in (("raw", g1), ("after p", g2), ("after pT", g3)):
        print(f"     {tag:<9s} V_l={g['V_l']:.15e}  m_l={g['m_l']:.15e}")

    # ---- 3. the volume pT_relax moves, and where it comes from --------
    print("=" * 78)
    print("3. volume exchanged (liquid receives, m^3)")
    dV_p = c_p.V_l - cell.V_l
    dV_pT = c_pT.V_l - cell.V_l
    print(f"   p_relax  alone            dV_l = {dV_p: .6e}")
    print(f"   pT_relax total            dV_l = {dV_pT: .6e}")
    print(f"   difference (thermal part) dV_l = {dV_pT - dV_p: .6e}"
          f"   ({abs((dV_pT-dV_p)/dV_pT)*100:.1f}% of the total)")
    print("   -> pT_relax computes the whole thing in one solve; nothing deferred.")

    # ---- 4. p_I closure insensitivity ---------------------------------
    print("=" * 78)
    print("4. p_relax closure sensitivity (p_I model)")
    base = None
    for mdl in ("gas", "liquid", "mean"):
        c = p_relax(cell, pI_model=mdl)
        if base is None:
            base = c.V_l - cell.V_l
        d = c.V_l - cell.V_l
        print(f"   p_I = {mdl:<7s} dV_l = {d: .10e}   rel. spread = "
              f"{abs(d/base-1):.3e}")

    # ---- 5. entropy increases along the sequence ----------------------
    print("=" * 78)
    print("5. entropy (J/K), must be non-decreasing")
    print(f"   initial   {r0['S']:.10f}")
    print(f"   p_relax   {c_p.report()['S']:.10f}")
    print(f"   pT_relax  {rT['S']:.10f}")
    print(f"   pTg_relax {rG['S']:.10f}")

    # pT and pTg are entropy maxima over their free coordinates
    # S(V,E,m) and g_Tp must share one entropy reference
    Vl,El,ml = cell.V_l, cell.E_l, cell.m_l
    Tl_, pl_ = liq.T(Vl,El,ml), liq.p(Vl,El,ml)
    s_ext = liq.S(Vl,El,ml)/ml
    s_gtp = (liq.cp*Tl_ + liq.b*pl_ + liq.q - liq.g_Tp(Tl_,pl_))/Tl_
    print(f"   liquid entropy reference consistent: "
          f"|s_ext - s_from_g| = {abs(s_ext-s_gtp):.3e} J/(kg K)")

    print("   pT_relax is a maximum over (V_l, E_l): sampling neighbours")
    worse = True
    for dv in (-1e-7, 1e-7):
        for de in (-1e-1, 1e-1):
            s = c_pT.report(c_pT.V_l + dv * cell.V, c_pT.E_l + de)['S']
            worse &= (s <= rT['S'])
    print(f"   all four neighbours have lower entropy: {worse}")

    # reduced entropy along m_l, each point at (p,T) equilibrium
    print("   Sigma(m_l) with (V_l,E_l) at pT equilibrium -- max must sit at pTg")
    base = c_pTg.m_l
    rows = []
    for f in (-4e-3,-2e-3,-1e-3,0.0,1e-3,2e-3,4e-3):
        ml_try = base*(1.0+f)
        c = cell.copy_with(cell.V_l, cell.E_l, ml_try)
        try:
            cc = pT_relax(c)
            rows.append((f, cc.report()['S'], cc.report()['F_C']))
        except Exception as ex:
            rows.append((f, float('nan'), float('nan')))
    smax = max(r[1] for r in rows)
    for f,S,FC in rows:
        mark = "  <== max" if S==smax else ""
        print(f"     dm_l/m_l={f:+.1e}   S={S:.8f}   F_C={FC: .4e}{mark}")
    print("=" * 78)


if __name__ == "__main__":
    main()
