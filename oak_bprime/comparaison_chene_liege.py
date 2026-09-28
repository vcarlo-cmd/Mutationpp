#!/usr/bin/env python
"""
Comparaison chêne / liège seul / liège-phénolique P50 à partir des tables B'
versionnées.

Aucun nouveau calcul d'équilibre : le script relit les CSV déjà produits par
oak_bprime.py, ../cork_bprime/cork_bprime.py (P50) et
../cork_bprime/cork_seul_bprime.py (liège sans résine).

    - chêne                : rendement en char 20 %, k = (1-y)/y = 4.0
    - liège seul           : rendement en char 20 % (imposé), k = 4.0
    - liège/phénolique P50 : rendement en char 20 % (TGA), k = 4.0

Les trois matériaux ont donc le même rendement en char, le même char (C pur)
et le même bord (air) : seule la composition du gaz de pyrolyse diffère.

Sorties :
    - tableaux imprimés (1 atm), repris dans comparaison_chene_liege.md
    - comparaison_chene_liege.csv : point stationnaire des trois matériaux
    - comparaison_chene_liege.png : B'c(T) à B'g fixé et point stationnaire
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CORK = os.path.join(HERE, "..", "cork_bprime")

Y_CHAR = 0.20                      # rendement en char, les trois matériaux
K = (1.0 - Y_CHAR) / Y_CHAR        # B'g/B'c stationnaire = 4.0

MATERIALS = {
    "chêne": {
        "dir": HERE, "prefix": "oak_bprime",
        "rho_v": 700.0,            # ordre de grandeur (oak_bprime/README.md)
        "gas": (0.221, 0.536, 0.243),
        "color": "tab:brown", "tag": "oak",
    },
    "liège seul": {
        "dir": CORK, "prefix": "cork_seul_bprime",
        "rho_v": 0.8 * 465.6,      # HYPOTHÈSE : P50 sans sa résine, volume constant
        "gas": (0.2576, 0.6133, 0.1291),
        "color": "tab:orange", "tag": "corkpure",
    },
    "liège P50": {
        "dir": CORK, "prefix": "cork_bprime",
        "rho_v": 465.6,            # mesurée (P50)
        "gas": (0.287, 0.592, 0.121),
        "color": "tab:green", "tag": "corkP50",
    },
}

P_BAR_1ATM = 1.01325
T_REPORT = (1000, 1500, 2000, 2500, 3000)
BG_REPORT = (0.0, 0.5, 2.0, 5.0)
BG_MAX_SWEEP = 10.0                # borne du balayage du point stationnaire


def read_csv(path):
    with open(path) as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def load(mat):
    d, p = mat["dir"], mat["prefix"]
    bc = read_csv(os.path.join(d, p + "_bc_table.csv"))
    ss = read_csv(os.path.join(d, p + "_steady_state.csv"))
    table = {}
    for r in bc:
        if abs(r["P_bar"] - P_BAR_1ATM) < 1e-6:
            table.setdefault(r["Bg"], {})[int(r["Tw_K"])] = r["Bc"]
    steady = {int(r["T_K"]): r for r in ss
              if abs(r["P_atm"] - 1.0) < 1e-9 and r["Bg_ss"] < BG_MAX_SWEEP}
    return table, steady


def main():
    names = list(MATERIALS)
    data = {name: load(m) for name, m in MATERIALS.items()}
    sep = " / "

    print("Gaz de pyrolyse (fractions molaires élémentaires)")
    for name, m in MATERIALS.items():
        c, h, o = m["gas"]
        print(f"  {name:10s} C {c:.3f}  H {h:.3f}  O {o:.3f}  O/C = {o / c:.2f}")
    print(f"\nRendement en char {Y_CHAR:.0%} pour les trois : k = {K:.1f}\n")

    print("B'c à 1 atm, B'g imposé (" + sep.join(names) + ")")
    for T in T_REPORT:
        cells = [sep.join(f"{data[n][0][bg][T]:.3f}" for n in names)
                 for bg in BG_REPORT]
        print(f"  {T:5d} | " + " | ".join(cells))

    print("\nPoint stationnaire (B'g = 4 B'c), 1 atm : B'c, B'g, "
          "(1+k)B'c/ρv [1e-3 m3/kg] (" + sep.join(names) + ")")
    rows = []
    for T in T_REPORT:
        ss = [data[n][1].get(T) for n in names]
        if any(s is None for s in ss):
            continue
        rec = [(1 + K) * s["Bc_ss"] / MATERIALS[n]["rho_v"]
               for n, s in zip(names, ss)]
        bc_ref = ss[-1]["Bc_ss"]
        print(f"  {T:5d} | B'c " + sep.join(f"{s['Bc_ss']:.4f}" for s in ss)
              + " | B'g " + sep.join(f"{s['Bg_ss']:.3f}" for s in ss)
              + " | rec " + sep.join(f"{1e3 * r:.3f}" for r in rec)
              + " | B'c/B'c(P50) " + sep.join(f"{s['Bc_ss'] / bc_ref:.2f}"
                                              for s in ss)
              + f" | B'c(B'g=0) {ss[0]['Bc_Bg0']:.4f}")
        rows.append([T, ss[0]["Bc_Bg0"]] + [s["Bc_ss"] for s in ss]
                    + [s["Bg_ss"] for s in ss] + [s["hw_ss_MJkg"] for s in ss]
                    + rec)

    out_csv = os.path.join(HERE, "comparaison_chene_liege.csv")
    tags = [MATERIALS[n]["tag"] for n in names]
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "Bc_Bg0"] + [f"Bc_ss_{t}" for t in tags]
                   + [f"Bg_ss_{t}" for t in tags]
                   + [f"hw_ss_{t}_MJkg" for t in tags]
                   + [f"recession_over_mdote_{t}_m3kg" for t in tags])
        for r in rows:
            w.writerow([r[0]] + [f"{v:.6e}" for v in r[1:]])
    print(f"\n-> {os.path.relpath(out_csv, HERE)}")

    # --- figure ------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    styles = {0.5: "--", 2.0: "-.", 5.0: "-"}
    t0 = sorted(data["chêne"][0][0.0])
    ax1.plot(t0, [data["chêne"][0][0.0][t] for t in t0], ":", color="k",
             label="B'g = 0 (les trois)")
    for name, m in MATERIALS.items():
        table, steady = data[name]
        for bg, ls in styles.items():
            Ts = sorted(table[bg])
            ax1.plot(Ts, [table[bg][t] for t in Ts], ls, color=m["color"],
                     label=f"{name}, B'g = {bg:g}")
        Ts = sorted(steady)
        ax2.plot(Ts, [steady[t]["Bc_ss"] for t in Ts], "-", color=m["color"],
                 label=f"{name} (k = {K:.0f})")
    Ts = sorted(data["chêne"][1])
    ax2.plot(Ts, [data["chêne"][1][t]["Bc_Bg0"] for t in Ts], ":", color="k",
             label="sans pyrolyse (B'g = 0)")
    for ax in (ax1, ax2):
        ax.set_xlim(300, 3500)
        ax.set_xlabel("T paroi [K]")
        ax.set_ylabel("B'c")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax1.set_ylim(0, 1.2)
    ax2.set_ylim(0, 0.8)
    ax1.set_title("B'c à B'g imposé, air, 1 atm")
    ax2.set_title("Point stationnaire B'g = 4 B'c (rendement char 20 %), 1 atm")
    fig.suptitle("Chêne / liège seul / liège-phénolique P50 — "
                 "même char, même k, gaz différents")
    fig.tight_layout()
    out_png = os.path.join(HERE, "comparaison_chene_liege.png")
    fig.savefig(out_png, dpi=150)
    print(f"-> {os.path.relpath(out_png, HERE)}")


if __name__ == "__main__":
    main()
