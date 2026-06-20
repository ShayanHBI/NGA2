#!/usr/bin/env python3
"""
JOINT NASG water fit (liquid + vapor) to CoolProp (IAPWS-95) reference data, over named regimes,
with cross-evaluation so the pre-impact-vs-broad compromise is explicit.

Targets (one consistent parameter set for mechanics AND phase change):
  mechanical : rho(T,P), c(T,P)           (CoolProp PropsSI over isobars)
  saturation : p_sat(T), latent heat L(T), rho_v(T), vapor cp   (CoolProp, Q=0/1)

gamma/cv are model knobs (not "real Cp/Cv"); only the vapor gas constant R_V is
physically pinned (= R_u/M_water). Gauge q_L=qp_L=0.

Deps: numpy, scipy, matplotlib. Usage: python3 fit_nasg.py
"""
import math, sys
import numpy as np
from scipy.optimize import least_squares
try:
    import CoolProp.CoolProp as CP
except ImportError:
    sys.exit("fit_nasg.py needs CoolProp:  python3 -m pip install CoolProp")

# ============================ CONFIG ============================
GammaG, Ms  = 1.3, 5.12
T1_amb, p1_amb, rho_air_amb = 288.15, 101325.0, 1.225   # ambient air

# Regimes to fit and compare. State-space (SG runs): pre-impact TL 264-348 K,
# PL -15..22 MPa; impact compresses to ~GPa (IAPWS data only to ~1 GPa / 1273 K).
SCENARIOS = {
    "pre-impact": dict(
        TEMPS=[274., 290., 305., 320., 335., 348.], P_LOW=1., P_HIGH=300., P_INC=15.,
        SAT_TLOW=274., SAT_THIGH=360.,
        PLOT_P=[1., 50., 150., 300.], PLOT_TLOW=274., PLOT_THIGH=352., PLOT_TINC=4.),
    "wideP": dict(    # cold liquid, extended pressure (as far as liquid exists before ice)
        TEMPS=[274., 290., 305., 320., 335., 348.], P_LOW=1., P_HIGH=6000., P_INC=200.,
        SAT_TLOW=274., SAT_THIGH=360.,
        PLOT_P=[1., 2000., 4000., 6000.], PLOT_TLOW=274., PLOT_THIGH=352., PLOT_TINC=4.),
}

W_RHO, W_C = 1.0, 1.0      # balanced rho/c (knee of the NASG Pareto frontier)
W_PSAT, W_L, W_RHOV, W_CPV = 1.0, 1.0, 0.5, 1.0
CEILING_MIN = 2000.0       # min NASG co-volume ceiling [kg/m3]; caps b at 1/CEILING_MIN so the impact (peak ~1947) stays below it
FIT_PHASES = {"liquid", "supercritical"}
DEPLOY = "pre-impact"      # scenario whose nondim EOS block is emitted for the input file
R_VAP = 8.31446 / 0.0180153

# Model-vs-NIST comparison plots: isobars [MPa] over a T range [K] (Tlow, Thigh, Tinc-for-dots)
PLOT_ISOBARS_MPA = [0.1, 1., 10., 20., 50., 100., 200., 500., 1000.]
PLOT_ISOBARS_T   = (280., 1000., 20.)
# Water principal Hugoniot, Rice & Walsh 1957 (J. Chem. Phys. 26, 824): tabulated (up, us) [m/s].
# rho0 = 1/V0 with V0 = 1.0018 cc/g -> 998.2 kg/m3; P0 = 1 atm.  (P[GPa] = rho0*us*up cross-checks the table.)
RICE_WALSH = dict(
    up=1.0e3*np.array([0.952, 1.392, 1.411, 1.655, 1.798, 1.806, 1.829, 2.370, 2.385, 4.13, 4.60, 4.72, 4.72, 4.81]),
    us=1.0e3*np.array([3.354, 4.093, 4.126, 4.536, 4.757, 4.777, 4.811, 5.601, 5.626, 8.07, 8.45, 8.49, 8.58, 8.74]),
    rho0=998.2, P0=1.0e5)
# Shayan's NASG fit over [300-500] K (SI), liquid + vapor -- overlaid for comparison.
# SHAYAN = dict(gL=1.1469, pinfL=5.5618e8, bL=6.8430e-4, cvL=3737.6, cpL=4286.5, qL=-1.1784e6, qpL=0.0,
#               gV=1.5372, cvV=856.75, cpV=1317.0, qV=2.1757e6, qpV=1.7877e4)
SHAYAN = dict(gL=1.4820e+00, pinfL=8.1001e+08, bL=4.7563e-04, cvL=2.9112e+03, cpL=4.3145e+03, qL=-1.1895e+06, qpL=0.0,
              gV=1.5367e+00, cvV=8.5759e+02, cpV=1.3179e+03, qV=2.1755e+06, qpV=3.0708e+02)
# ===============================================================


def _phase_id(T, P):
    """CoolProp phase at (T,P) -> 0=liquid, 1=dense supercritical, 2=vapor/low-density gas, -1=two-phase.
    supercritical_gas (T>Tc, P<Pc) is dilute steam -> classed as vapor, not dense supercritical."""
    ph = CP.PhaseSI('T', T, 'P', P, 'Water')
    if ph in ('gas', 'vapor', 'supercritical_gas'):  return 2
    if 'supercritical' in ph:                        return 1   # supercritical / supercritical_liquid (dense)
    if ph == 'liquid':                               return 0
    return -1


def wagner_psat(T):
    """IAPWS-95 saturation pressure [Pa], valid 273.16-647.096 K."""
    Tc, pc = 647.096, 22.064e6
    a = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    e = (1, 1.5, 3, 3.5, 4, 7.5)
    tau = 1.0 - T / Tc
    s = sum(a[i] * tau ** e[i] for i in range(6))
    return pc * np.exp(Tc / T * s)


def gather_liquid(scn):
    cols = {k: [] for k in ("P", "T", "rho", "c")}
    Ps = np.arange(scn["P_LOW"], scn["P_HIGH"] + 1e-9, scn["P_INC"]) * 1e5   # bar -> Pa
    for Ti in scn["TEMPS"]:
        kept = 0
        for P in Ps:
            if _phase_id(Ti, P) not in (0, 1):   # liquid / supercritical only (FIT_PHASES)
                continue
            cols["T"].append(Ti); cols["P"].append(P)
            cols["rho"].append(CP.PropsSI('D', 'T', Ti, 'P', P, 'Water'))
            cols["c"].append(CP.PropsSI('A', 'T', Ti, 'P', P, 'Water')); kept += 1
        print(f"    liquid T={Ti:6.1f} K : {kept} pts")
    return {k: np.array(v) for k, v in cols.items()}


def gather_sat(scn):
    psat = lambda T: CP.PropsSI('P', 'T', T, 'Q', 0, 'Water')          # saturation pressure
    psT = np.linspace(scn["SAT_TLOW"], min(scn["SAT_THIGH"], 646.0), 35)
    wT  = np.linspace(scn["SAT_TLOW"], min(scn["SAT_THIGH"], 640.0), 25)
    sat = lambda f, Q, T: np.array([CP.PropsSI(f, 'T', t, 'Q', Q, 'Water') for t in T])
    print(f"    saturation : {len(psT)} p_sat + {len(wT)} L/rhov/cpv  ({psT[0]:.0f}-{psT[-1]:.0f} K)  [CoolProp]")
    return dict(psT=psT, psat=np.array([psat(t) for t in psT]),
                wT=wT, psw=np.array([psat(t) for t in wT]),
                rho_v=sat('D', 1, wT), L=sat('H', 1, wT) - sat('H', 0, wT), cpv=sat('C', 1, wT))


def fetch_isobar(scn, P_bar):
    a = fetch_isobar_phases(scn, P_bar)
    return a[np.isin(a[:, 3], (0, 1))][:, :3] if len(a) else a   # liquid / supercritical only


def fetch_isobar_phases(scn, P_bar):
    """Isobar over the plot T-range; returns (T, rho, c, phase_id): 0=liquid, 1=supercritical, 2=vapor."""
    P = P_bar * 1e5
    out = []
    for T in np.arange(scn["PLOT_TLOW"], scn["PLOT_THIGH"] + 1e-9, scn["PLOT_TINC"]):
        pid = _phase_id(T, P)
        if pid < 0:
            continue
        out.append((T, CP.PropsSI('D', 'T', T, 'P', P, 'Water'), CP.PropsSI('A', 'T', T, 'P', P, 'Water'), float(pid)))
    return np.array(out)


# --- model ---
def rho_nasg(P, T, g, pinf, b, cv): return (P + pinf) / ((g - 1.) * cv * T + b * (P + pinf))
def c_nasg(P, T, g, pinf, b, cv):
    r = rho_nasg(P, T, g, pinf, b, cv); return np.sqrt(np.maximum(0., g * (P + pinf) / (r * (1. - b * r))))
def e_nasg(P, rho, g, pinf, b, cv):  # internal energy along the liquid EOS (q_L=0 gauge)
    v = 1.0/rho; return (P + g*pinf)*(v - b)/(g - 1.0)
def hugoniot_nasg(d, rho, P0, rho0):
    """NASG principal Hugoniot P(rho): Rankine-Hugoniot e-e0 = 0.5(P+P0)(v0-v), closed form."""
    g, pinf, b = d['gL'], d['pinfL'], d['bL']
    v, v0 = 1.0/rho, 1.0/rho0
    e0  = e_nasg(P0, rho0, g, pinf, b, d['cvL'])
    num = e0 + 0.5*P0*(v0 - v) - g*pinf*(v - b)/(g - 1.0)
    den = (v - b)/(g - 1.0) - 0.5*(v0 - v)
    return num/den
def hugoniot_marsh(up, c0, s, rho0, P0):
    """us-up linear Hugoniot (Marsh) -> (rho, P) for a shock from (rho0, P0, at rest)."""
    us = c0 + s*up; return rho0*us/(us - up), P0 + rho0*us*up

def coef(p):
    gL, piL, bL, cvL, cvV, qV, qpV = p
    cpV, cpL = cvV + R_VAP, gL * cvL
    return (cpL-cpV+qpV)/R_VAP, -qV/R_VAP, (cpV-cpL)/R_VAP, (cpL-cvL)/R_VAP, bL/R_VAP, piL, cpL, cpV

def psat_resid(p, T, ps):
    AS, BS, CS, DS, ES, piL, *_ = coef(p)
    return AS + (BS + ES*ps)/T + CS*np.log(T) + DS*np.log(ps+piL) - np.log(ps)
def psat_curve(p, T):
    ps = wagner_psat(T).copy()
    for _ in range(80):
        ps = np.exp(psat_resid(p, T, ps) + np.log(ps))
    return ps
def L_model(p, T, ps):
    *_, cpL, cpV = coef(p); return (cpV - cpL)*T - p[2]*ps + p[5]
def rhov_model(T, ps): return ps / (R_VAP * T)

def L_cc(T):
    """Clausius-Clapeyron latent heat from Wagner p_sat (ideal vapor, v_l neglected;
    accurate cold where vapor is near-ideal, degrades approaching critical)."""
    dT = 0.02
    dpdT = (wagner_psat(T + dT) - wagner_psat(T - dT)) / (2 * dT)
    return T * (R_VAP * T / wagner_psat(T)) * dpdT

def shayan_sat(T):
    """Saturation p_sat(T), latent heat L(T), vapor density rho_v(T) for Shayan's NASG fit,
    using his own gauge (q_L != 0) and his vapor R_V = cpV - cvV."""
    s = SHAYAN
    RV = s['cpV'] - s['cvV']
    AS = (s['cpL'] - s['cpV'] + s['qpV'] - s['qpL'])/RV
    BS = (s['qL'] - s['qV'])/RV
    CS = (s['cpV'] - s['cpL'])/RV
    DS = (s['cpL'] - s['cvL'])/RV
    ES = s['bL']/RV
    ps = wagner_psat(T).copy()
    for _ in range(80):
        ps = np.exp(AS + (BS + ES*ps)/T + CS*np.log(T) + DS*np.log(ps + s['pinfL']))
    L = (s['cpV'] - s['cpL'])*T - s['bL']*ps + (s['qV'] - s['qL'])
    return ps, L, ps/(RV*T)


def pack(d): return [d['gL'], d['pinfL'], d['bL'], d['cvL'], d['cvV'], d['qV'], d['qpV']]

def evaluate(p, Ld, Sd):
    rms = lambda r: 100*float(np.sqrt(np.mean(np.asarray(r)**2)))
    g, pi, b, cv = p[0], p[1], p[2], p[3]
    return dict(
        rho=rms((rho_nasg(Ld['P'], Ld['T'], g, pi, b, cv)-Ld['rho'])/Ld['rho']),
        c  =rms((c_nasg(Ld['P'], Ld['T'], g, pi, b, cv)-Ld['c'])/Ld['c']),
        psat=rms(np.expm1(psat_resid(p, Sd['psT'], Sd['psat']))),
        L  =rms((L_model(p, Sd['wT'], Sd['psw'])-Sd['L'])/Sd['L']),
        rhov=rms((rhov_model(Sd['wT'], Sd['psw'])-Sd['rho_v'])/Sd['rho_v']))


def fit(Ld, Sd, W=None):
    if W is None: W = (W_RHO, W_C, W_PSAT, W_L, W_RHOV, W_CPV)
    wr, wc, wp, wl, wv, wcp = W
    def block(r, w): return w*np.asarray(r)/math.sqrt(max(1, len(r)))
    def res(p):
        g, pi, b, cv = p[0], p[1], p[2], p[3]
        rr = (rho_nasg(Ld['P'], Ld['T'], g, pi, b, cv)-Ld['rho'])/Ld['rho']
        rc = (c_nasg(Ld['P'], Ld['T'], g, pi, b, cv)-Ld['c'])/Ld['c']
        rp = psat_resid(p, Sd['psT'], Sd['psat'])
        rl = (L_model(p, Sd['wT'], Sd['psw'])-Sd['L'])/Sd['L']
        rv = (rhov_model(Sd['wT'], Sd['psw'])-Sd['rho_v'])/Sd['rho_v']
        rcp = ((p[4]+R_VAP)-Sd['cpv'])/Sd['cpv']
        return np.concatenate([block(rr, wr), block(rc, wc), block(rp, wp),
                               block(rl, wl), block(rv, wv), block(rcp, wcp)])
    x0 = [1.4, 8e8, 4e-4, 3500., 1500., 2.5e6, 1e4]
    lo = [1.01, 0., 0., 100., 200., 0., -5e4]
    hi = [8., 5e10, 1.0/CEILING_MIN, 1e4, 3000., 6e6, 5e4]   # b <= 1/CEILING_MIN keeps 1/b above the impact density
    p = least_squares(res, x0, bounds=(lo, hi), xtol=1e-14, ftol=1e-14, max_nfev=30000).x
    return dict(gL=p[0], pinfL=p[1], bL=p[2], cvL=p[3], cvV=p[4],
                gV=(p[4]+R_VAP)/p[4], qV=p[5], qpV=p[6])


def ref_scales():
    rhoG1 = ((GammaG-1)*Ms**2+2)/((GammaG+1)*Ms**2)
    pG1n  = 0.25*rhoG1/GammaG*((GammaG+1)*Ms/(Ms**2-1))**2
    CvGn  = (pG1n*(2*GammaG/(GammaG+1)*(Ms**2-1)+1))/(GammaG-1)
    TG1n  = pG1n/(CvGn*rhoG1*(GammaG-1))
    p_ref, rho_ref = p1_amb/pG1n, rho_air_amb/rhoG1
    return rho_ref, math.sqrt(p_ref/rho_ref), T1_amb/TG1n, p_ref


def nondim(d, sc):
    rho_ref, u_ref, T_ref, p_ref = sc; sT = T_ref/u_ref**2
    shift = lambda cp, R: cp*sT*math.log(T_ref) - R*sT*math.log(p_ref)
    RL = (d['gL']-1)*d['cvL']
    return dict(GammaL=d['gL'], PinfL=d['pinfL']/p_ref, b=d['bL']*rho_ref, cvL=d['cvL']*sT,
                qpL=shift(d['gL']*d['cvL'], RL), GammaV=d['gV'], cvV=d['cvV']*sT,
                qV=d['qV']/u_ref**2, qpV=d['qpV']*sT + shift(d['cvV']+R_VAP, R_VAP))


def emit_input(name, d, nd, sc, e):
    """Paste-ready nondimensional EOS block for amrcomp_impact/input, with a full
    summary of the assumptions that set the nondimensionalization. Keys match what
    simulation.f90 reads, so it can be pasted directly into the input EOS block."""
    rho_ref, u_ref, T_ref, p_ref = sc
    kv  = lambda k, v: f"{k:<18s}{v}"
    bar = "# " + "=" * 68
    return "\n".join([
        bar,
        "# NASG water EOS (nondimensional) -- paste into the EOS block of amrcomp_impact/input",
        f"# generated by fit_nasg.py   scenario: {name}   fit phases: {sorted(FIT_PHASES)}",
        "# -- assumptions feeding the nondimensionalization --",
        f"#   shock      Ms = {Ms} ,  GammaG = {GammaG}",
        f"#   ambient    p1 = {p1_amb:.6g} Pa ,  rho1 = {rho_air_amb:.6g} kg/m3 ,  T1 = {T1_amb:.6g} K",
        f"#   ref scales rho_ref = {rho_ref:.6g} kg/m3 ,  u_ref = {u_ref:.6g} m/s",
        f"#              T_ref   = {T_ref:.6g} K ,       p_ref = {p_ref:.6g} Pa",
        f"#   vapor R_V  {R_VAP:.6g} J/kg/K  (pinned = R_u / M_water)",
        "#   gauge      liquid q = qp = 0 in SI;  nondim qpL absorbs +cp*ln(T_ref) - R*ln(p_ref)",
        f"#   fit RMS    rho {e['rho']:.2f}% ,  c {e['c']:.2f}% ,  p_sat {e['psat']:.2f}% ,  L {e['L']:.2f}%",
        f"#   co-volume  rho_max = 1/b = {1.0/d['bL']:.0f} kg/m3  (rho_hat {1.0/d['bL']/rho_ref:.0f})",
        "# -- air (ideal gas: assumed, NOT fit -- must match the shock setup above) --",
        kv("GammaG:", f"{GammaG:.8g}"),
        "# -- liquid water (NASG); liquid q omitted = 0 by gauge (hardcoded in sim) --",
        kv("GammaL:",          f"{nd['GammaL']:.8g}"),
        kv("Liquid pinf:",     f"{nd['PinfL']:.8g}"),
        kv("Liquid covolume:", f"{nd['b']:.8g}"),
        kv("Liquid cv:",       f"{nd['cvL']:.8g}"),
        kv("Liquid qp:",       f"{nd['qpL']:.8g}"),
        "# -- water vapor (ideal gas) --",
        kv("GammaV:",          f"{nd['GammaV']:.8g}"),
        kv("Vapor cv:",        f"{nd['cvV']:.8g}"),
        kv("Vapor q:",         f"{nd['qV']:.8g}"),
        kv("Vapor qp:",        f"{nd['qpV']:.8g}"),
        bar,
    ])


def plot_mech(scn, name, d, fn):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    iso = {Pb: fetch_isobar(scn, Pb) for Pb in scn["PLOT_P"]}
    iso = {Pb: a for Pb, a in iso.items() if len(a)}
    Ps = sorted(iso); nP = len(Ps)
    fig, ax = plt.subplots(nP, 2, figsize=(11, 2.4*nP), sharex=True, squeeze=False)
    for row, Pb in enumerate(Ps):
        a = iso[Pb]; T = a[:, 0]; P = Pb*1e5
        for col, (yd, ylab, fnk) in enumerate([(a[:, 1], r"$\rho$ [kg/m³]", rho_nasg),
                                               (a[:, 2], "c [m/s]", c_nasg)]):
            ax_ = ax[row, col]
            yn = fnk(P, T, d['gL'], d['pinfL'], d['bL'], d['cvL'])
            ax_.plot(T, yd, 'k.', ms=5, label="CoolProp", zorder=3); ax_.plot(T, yn, 'C0-', lw=1.8, label="NASG")
            ax_.set_ylabel(ylab)
            ax_.set_title(f"{Pb:.0f} bar  (RMS {100*np.sqrt(np.mean(((yn-yd)/yd)**2)):.2f}%)", fontsize=9)
            if row == 0 and col == 0: ax_.legend(fontsize=8)
    for c in range(2): ax[-1, c].set_xlabel("T [K]")
    fig.suptitle(f"Mechanical [{name}]: rho, c vs T per isobar", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.99]); fig.savefig(fn, dpi=120); print(f"    {fn}")


def plot_sat(name, Sd, d, fn):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    p = pack(d)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    psS, LS, rhovS = shayan_sat(Sd['psT'])
    ax[0].semilogy(Sd['psT'], Sd['psat']/1e5, 'k.', ms=5, label="CoolProp")
    ax[0].semilogy(Sd['psT'], psat_curve(p, Sd['psT'])/1e5, 'C0-', label="NASG (mine)")
    ax[0].semilogy(Sd['psT'], psS/1e5, 'C1-.', label="NASG (Shayan)")
    ax[0].set(xlabel="T [K]", ylabel="p_sat [bar]", title="Saturation pressure"); ax[0].legend()
    ax[1].plot(Sd['psT'], L_model(p, Sd['psT'], Sd['psat'])/1e3, 'C0-', label="NASG (mine)")
    ax[1].plot(Sd['psT'], LS/1e3, 'C1-.', label="NASG (Shayan)")
    ax[1].plot(Sd['psT'], L_cc(Sd['psT'])/1e3, 'C2:', lw=2, label="Clausius-Clapeyron (IAPWS)")
    ax[1].plot(Sd['wT'], Sd['L']/1e3, 'k.', ms=5, label="CoolProp")
    ax[1].set(xlabel="T [K]", ylabel="L [kJ/kg]", title="Latent heat"); ax[1].legend()
    ax[2].semilogy(Sd['psT'], rhov_model(Sd['psT'], Sd['psat']), 'C0-', label="NASG (mine)")
    ax[2].semilogy(Sd['psT'], rhovS, 'C1-.', label="NASG (Shayan)")
    ax[2].semilogy(Sd['wT'], Sd['rho_v'], 'k.', ms=5, label="CoolProp")
    ax[2].set(xlabel="T [K]", ylabel=r"$\rho_v$ [kg/m³]", title="Sat. vapor density"); ax[2].legend()
    fig.suptitle(f"Saturation [{name}]: NASG vs CoolProp (IAPWS-95)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(fn, dpi=120); print(f"    {fn}")


def plot_isobars(d, fn):
    """Model-vs-NIST: c(T), rho(T) along isobars. NASG lines; NIST circles (liquid/supercritical), triangles (vapor)."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from matplotlib.lines import Line2D
    Tlo, Thi, Tinc = PLOT_ISOBARS_T
    Tline = np.linspace(Tlo, Thi, 300)
    norm = colors.LogNorm(vmin=min(PLOT_ISOBARS_MPA), vmax=max(PLOT_ISOBARS_MPA)); cmap = cm.viridis
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    for Pm in PLOT_ISOBARS_MPA:
        col = cmap(norm(Pm)); Ppa = Pm*1e6
        ax[0].plot(Tline, c_nasg  (Ppa, Tline, d['gL'], d['pinfL'], d['bL'], d['cvL']), '-', color=col, lw=1.6)
        ax[1].plot(Tline, rho_nasg(Ppa, Tline, d['gL'], d['pinfL'], d['bL'], d['cvL']), '-', color=col, lw=1.6)
        s = SHAYAN
        ax[0].plot(Tline, c_nasg  (Ppa, Tline, s['gL'], s['pinfL'], s['bL'], s['cvL']), '-.', color=col, lw=1.3)
        ax[1].plot(Tline, rho_nasg(Ppa, Tline, s['gL'], s['pinfL'], s['bL'], s['cvL']), '-.', color=col, lw=1.3)
        a = fetch_isobar_phases({"PLOT_TLOW": Tlo, "PLOT_THIGH": Thi, "PLOT_TINC": Tinc}, Pm*10.0)  # MPa -> bar
        if not len(a):
            continue
        for pid, mk in ((0, ''), (1, 'o'), (2, '^')):   # NIST: liquid=line only, supercritical=circle, vapor=triangle
            m = a[:, 3] == pid
            if m.any():
                ax[0].plot(a[m, 0], a[m, 2], ls='--', marker=mk, color=col, lw=1.0, ms=3.6, mec='k', mew=0.3, zorder=3)
                ax[1].plot(a[m, 0], a[m, 1], ls='--', marker=mk, color=col, lw=1.0, ms=3.6, mec='k', mew=0.3, zorder=3)
        # Connect liquid <-> vapor across the saturation jump: vertical line at ~T_sat (sub-critical isobars)
        liq, vap = a[a[:, 3] == 0], a[a[:, 3] == 2]
        if len(liq) and len(vap):
            lp, vp = liq[np.argmax(liq[:, 0])], vap[np.argmin(vap[:, 0])]   # last liquid, first vapor
            Tj = 0.5*(lp[0] + vp[0])
            ax[0].plot([Tj, Tj], [lp[2], vp[2]], ls=':', color=col, lw=1.0, zorder=2)
            ax[1].plot([Tj, Tj], [lp[1], vp[1]], ls=':', color=col, lw=1.0, zorder=2)
    ax[0].set(xlabel="T [K]", ylabel="c [m/s]",         title="Sound speed")
    ax[1].set(xlabel="T [K]", ylabel=r"$\rho$ [kg/m³]", title="Density")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30); cb.set_label("P [MPa]")
    ax[0].legend(handles=[Line2D([], [], color='k', ls='-',  label='NASG (mine)'),
                          Line2D([], [], color='k', ls='-.', label='NASG (Shayan)'),
                          Line2D([], [], color='k', ls='--', label='CoolProp liquid'),
                          Line2D([], [], color='k', ls='--', marker='o', ms=4, label='CoolProp supercritical'),
                          Line2D([], [], color='k', ls='--', marker='^', ms=4, label='CoolProp vapor')], fontsize=7)
    fig.suptitle("NASG: mine (solid) vs Shayan (dash-dot) vs CoolProp (dashed+markers), isobars 0.1–1000 MPa", fontsize=10)
    fig.savefig(fn, dpi=120, bbox_inches='tight'); print(f"    {fn}")


def plot_hugoniot(d, fn):
    """Principal Hugoniot vs piston (particle) velocity: Ms, P2 [GPa], rho2/rho1 = f(u_p).
    NASG (closed form) vs Rice-Walsh 1957 water data."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    g, pinf, b = d['gL'], d['pinfL'], d['bL']
    rho0, P0 = RICE_WALSH['rho0'], RICE_WALSH['P0']
    a0 = np.sqrt(g*(P0+pinf)/(rho0*(1.-b*rho0)))            # ambient bulk sound speed [m/s]
    # NASG locus, parametrized by compressed density up to the co-volume ceiling
    def locus(dd):   # NASG principal Hugoniot -> (up, Ms, P2[GPa], rho2/rho1), common a0
        r2 = np.linspace(rho0*1.0005, 0.98/dd['bL'], 400)
        P2_ = hugoniot_nasg(dd, r2, P0, rho0)
        up_ = np.sqrt(np.maximum(0., (P2_-P0)*(1./rho0 - 1./r2)))
        us_ = (P2_-P0)/(rho0*np.maximum(up_, 1e-12))
        k = (P2_ > P0) & (up_ > 0)
        return up_[k], us_[k]/a0, P2_[k]/1e9, r2[k]/rho0
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
    ylims = [(1., 5.), (0., 12.), (1.0, 2.1)]   # Ms, P2 [GPa], rho2/rho1
    for lab, dd, sty in [("NASG (mine)", d, dict(ls='-', color='C0', lw=1.9)),
                         ("NASG (Shayan)", SHAYAN, dict(ls='-.', color='C1', lw=1.6))]:
        u_, Ms_, P2g_, rr_ = locus(dd)
        for axi, y in zip(ax, [Ms_, P2g_, rr_]):
            axi.plot(u_, y, label=lab, **sty)
    for axi, ylab, yl in zip(ax, [r"$M_s = u_s/a_0$", r"$P_2$ [GPa]", r"$\rho_2/\rho_1$"], ylims):
        axi.set(xlabel=r"$u_\mathrm{piston}$ [m/s]", ylabel=ylab, xlim=(0., 3000.))
        if yl is not None:
            axi.set_ylim(*yl)
    upM, usM = RICE_WALSH['up'], RICE_WALSH['us']
    for axi, y in zip(ax, [usM/a0, (P0+rho0*usM*upM)/1e9, usM/(usM-upM)]):
        axi.plot(upM, y, 'ks', ms=5, mfc='none', label="Rice-Walsh 1957")
    ax[0].legend(fontsize=8)
    fig.suptitle(f"Principal Hugoniot vs piston velocity  (a0 = {a0:.0f} m/s)", fontsize=11)
    fig.savefig(fn, dpi=120, bbox_inches='tight'); print(f"    {fn}")


def main():
    fits, data = {}, {}
    for name, scn in SCENARIOS.items():
        print(f"Fetching + fitting [{name}] ...")
        Ld, Sd = gather_liquid(scn), gather_sat(scn)
        data[name] = (Ld, Sd); fits[name] = fit(Ld, Sd)

    rho_ref = ref_scales()[0]
    print("\n=== Fitted NASG parameters (SI) ===")
    for name, d in fits.items():
        print(f"  [{name:10s}] liq: g={d['gL']:.4f} pinf={d['pinfL']:.3e} b={d['bL']:.3e} cv={d['cvL']:.0f}"
              f"  vap: g={d['gV']:.4f} cv={d['cvV']:.0f} qV={d['qV']:.3e} qpV={d['qpV']:.3e}")
    print("\n=== Co-volume density ceiling  (impact reaches rho ~1722 kg/m3 = rho_hat ~230) ===")
    for name, d in fits.items():
        rmax = 1.0 / d['bL']
        print(f"  [{name:10s}] rho_max = {rmax:7.0f} kg/m3  (rho_hat {rmax/rho_ref:5.0f})  -> impact at {1722/rmax*100:.0f}% of ceiling")

    print("\n=== Cross-evaluation: RMS%% of each fit on each regime's data ===")
    print("  fit \\ data    " + "".join(f"{n:>26s}" for n in SCENARIOS))
    for fn_, d in fits.items():
        row = f"  {fn_:14s}"
        for dn, (Ld, Sd) in data.items():
            e = evaluate(pack(d), Ld, Sd)
            row += f"  rho{e['rho']:4.1f} c{e['c']:4.1f} ps{e['psat']:4.1f} L{e['L']:4.1f}"
        print(row)

    # rho/c trade: weight sweep on the pre-impact regime
    print("\n=== Pre-impact rho/c trade (weight sweep) ===")
    Ld, Sd = data["pre-impact"]
    #          (W_rho,W_c,W_psat,W_L,W_rhov,W_cpv)
    presets = {
        "c-favored  (Wc=3)":  (1, 3, 1, 1, .5, 1),
        "balanced   (Wc=1)":  (1, 1, 1, 1, .5, 1),
        "rho-favored(Wr=3)":  (3, 1, 1, 1, .5, 1),
        "rho-heavy  (Wr=10)": (10, 1, 1, 1, .5, 1),
        "rho-only   (Wr=30)": (30, 1, 1, 1, .5, 1),
        "sat-light":          (1, 1, .2, .2, .1, .2),
        "mech-only (no sat)": (1, 1, 0, 0, 0, 0),
    }
    print(f"  (co-volume ceiling capped at 1/b >= {CEILING_MIN:.0f} kg/m3)")
    print(f"  {'preset':20s} {'rho%':>6s} {'c%':>6s} {'psat%':>7s} {'L%':>6s}   gL     b         1/b")
    for name, W in presets.items():
        d = fit(Ld, Sd, W); e = evaluate(pack(d), Ld, Sd)
        print(f"  {name:20s} {e['rho']:6.2f} {e['c']:6.2f} {e['psat']:7.2f} {e['L']:6.2f}"
              f"   {d['gL']:.3f} {d['bL']:.2e} {1.0/d['bL']:5.0f}")

    sc = ref_scales()
    print(f"\n=== Reference scales ===  rho={sc[0]:.3f}  u={sc[1]:.1f}  T={sc[2]:.1f}  p={sc[3]:.3e}")
    print("\n=== Nondimensional coefficients (per scenario) ===")
    for name, d in fits.items():
        nd = nondim(d, sc)
        print(f"  [{name}]")
        print(f"    Liquid: GammaL {nd['GammaL']:.5f}  PinfL {nd['PinfL']:.4f}  b {nd['b']:.6f}  "
              f"cvL {nd['cvL']:.5f}  qL 0  qpL {nd['qpL']:.4f}")
        print(f"    Vapor : GammaV {nd['GammaV']:.5f}  cvV {nd['cvV']:.5f}  qV {nd['qV']:.5f}  qpV {nd['qpV']:.4f}")

    # Paste-ready input EOS block for the deployed scenario (single source of truth)
    d = fits[DEPLOY]; Ld, Sd = data[DEPLOY]
    print("\n" + emit_input(DEPLOY, d, nondim(d, sc), sc, evaluate(pack(d), Ld, Sd)))

    try:
        for name, scn in SCENARIOS.items():
            Ld, Sd = data[name]
            plot_mech(scn, name, fits[name], f"fit_nasg_mech_{name}.png")
            plot_sat(name, Sd, fits[name], f"fit_nasg_sat_{name}.png")
        plot_isobars(fits[DEPLOY], "fit_nasg_isobars.png")
        plot_hugoniot(fits[DEPLOY], "fit_nasg_hugoniot.png")
    except Exception as e:
        print(f"(plotting skipped: {e})")


if __name__ == "__main__":
    main()
