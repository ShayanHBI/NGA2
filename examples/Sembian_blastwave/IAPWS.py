import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# IAPWS / Wagner-Pruss saturation pressure equation for ordinary water
# T in K, p_sat in Pa

Tc = 647.096      # K
Pc = 22.064e6     # Pa

a = np.array([
    -7.85951783,
     1.84408259,
    -11.7866497,
     22.6807411,
    -15.9618719,
     1.80122502,
])

exps = np.array([1.0, 1.5, 3.0, 3.5, 4.0, 7.5])

def psat_iapws(T):
    T = np.asarray(T, dtype=float)
    theta = 1.0 - T / Tc

    series = np.zeros_like(T)
    for ai, ei in zip(a, exps):
        series += ai * theta**ei

    return Pc * np.exp((Tc / T) * series)

# Temperature range
T = np.linspace(300.0, 500.0, 50)
p = psat_iapws(T)

df = pd.DataFrame({
    "T_K": T,
    "T_C": T - 273.15,
    "p_sat_Pa": p,
    "p_sat_bar": p / 1e5,
    "p_sat_MPa": p / 1e6,
})

out_dir = Path(".")

csv_path = out_dir / "iapws_water_saturation_p_vs_T.csv"
png_path = out_dir / "iapws_water_saturation_p_vs_T.png"
pdf_path = out_dir / "iapws_water_saturation_p_vs_T.pdf"

df.to_csv(csv_path, index=False)

plt.figure(figsize=(7.0, 4.8))
plt.semilogy(df["T_K"], df["p_sat_Pa"], linewidth=2)
plt.xlabel("Temperature, T [K]")
plt.ylabel("Saturation pressure, $p_{sat}$ [Pa]")
plt.title("Pure-water liquid-vapor saturation curve")
plt.grid(True, which="both", alpha=0.35)
plt.tight_layout()
plt.savefig(png_path, dpi=220)
plt.savefig(pdf_path)

# Sparse checkpoint table
T_sparse = np.array([
    300, 320, 340, 360, 373.15,
    380, 400, 420, 450, 480, 500
])

p_sparse = psat_iapws(T_sparse)

sparse = pd.DataFrame({
    "T_K": T_sparse,
    "T_C": T_sparse - 273.15,
    "p_sat_Pa": p_sparse,
    "p_sat_bar": p_sparse / 1e5,
    "p_sat_MPa": p_sparse / 1e6,
})

sparse_path = out_dir / "iapws_water_saturation_p_vs_T_sparse.csv"
sparse.to_csv(sparse_path, index=False)

print(sparse)