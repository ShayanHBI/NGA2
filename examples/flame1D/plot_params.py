import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 22})
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["legend.handlelength"] = 0.8


lw = 2.5
ms = 8
mk = ["", "s", "x", "o", "^", "+", "*", "D", ">", "o", "o", "o", "o", "o", "o"]
ls = ["-", "--", "-.", ":"]

c = [
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    # "#ffff33",
    "#a65628",
    "#f781bf",
]
c = np.array(c)

