import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# File paths
ref_path   = Path("./iapws_water_saturation_p_vs_T.csv")
pure_path  = Path("./ptg_relax_pT_pure_water.csv")
mix_path   = Path("./ptg_relax_pT_Yv0_0p5.csv")

# Read data
ref  = pd.read_csv(ref_path)
pure = pd.read_csv(pure_path)
mix  = pd.read_csv(mix_path)

# Keep only valid states if that column exists
if "valid_state" in pure.columns:
    pure = pure[pure["valid_state"] == 1]

if "valid_state" in mix.columns:
    mix = mix[mix["valid_state"] == 1]

# Create plot
plt.figure(figsize=(7.5, 5.0))

# Reference saturation curve
plt.semilogy(
    ref["T_K"],
    ref["p_sat_Pa"],
    linewidth=2,
    label="IAPWS"
)

# PTg_relax pure-water results
plt.semilogy(
    pure["T_final"],
    pure["pv_final"],
    "o",
    markersize=6,
    label="Liquid-vapor (Pelanti and Shyue 2014 algorithm into NGA2)"
)

# PTg_relax mixture results
plt.semilogy(
    mix["T_final"],
    mix["pv_final"],
    "s",
    markersize=3,
    label="Liquid-vapor-air (NGA2)"
)

plt.xlabel("Temperature, T [K]")
plt.ylabel("Pressure [Pa]")
plt.title("PTg_relax results vs IAPWS reference saturation curve")
plt.grid(True, which="both", alpha=0.35)
plt.legend()
plt.tight_layout()

# Save
pdf_path = Path("./ptg_relax_vs_iapws_pT.pdf")

plt.savefig(pdf_path)

print("Saved:")
print(pdf_path)