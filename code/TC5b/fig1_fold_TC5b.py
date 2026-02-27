#! /usr/bin/env python3
"""
File:           fig1_fold_TC5b.py
Description:    Fold TC5b in the HP-model and visualize the conformation.
"""
import os
from pathlib import Path

import scienceplots  # noqa: F401

from matplotlib import pyplot as plt
from prospr import Protein, depth_first_bnb, plot_protein

os.chdir(Path(__file__).parent)

TC5b = Protein("PHHHPHHPPHHHPPHPHHHP", dim=2)

print("Folding TC5b in 2D:")
depth_first_bnb(TC5b)
moves = ["v", "<", None, ">", "^"]
print(f"Conformation: {[moves[i + 2] for i in TC5b.hash_fold()]}")
print(f"Stability: {TC5b.score}")

plt.style.use(["science"])
plt.rcParams.update({
    "legend.frameon": True,
    "legend.framealpha": 1,
    "legend.facecolor": "white",
    "legend.edgecolor": "black",
})
fig = plt.figure()
plot_protein(TC5b, style="paper", show=False, ax=plt.gca(), annotate_first=True)

plt.tight_layout()
fig = plt.gcf()
w, h = fig.get_size_inches()
fig.set_size_inches(w * 1.25, h * 1.3)
plt.savefig("fig1_TC5b.eps")
plt.savefig("fig1_TC5b.png")
