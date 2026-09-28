#!/usr/bin/env python
"""
Comparaison chêne / liège-phénolique P50 à partir des tables B' versionnées.

Aucun nouveau calcul d'équilibre : le script relit les CSV déjà produits par
oak_bprime.py et ../cork_bprime/cork_bprime.py.

    - chêne              : rendement en char 20 %, k = (1-y)/y = 4.0
    - liège/phénolique P50 : rendement en char 20 % (TGA), k = 4.0

Les deux matériaux ont donc le même rendement en char, le même char (C pur)
et le même bord (air) : seule la composition du gaz de pyrolyse diffère.

Sorties :
    - tableaux imprimés (1 atm), repris dans comparaison_chene_liege.md
    - comparaison_chene_liege.csv : point stationnaire des deux matériaux
    - comparaison_chene_liege.png : B'c(T) à B'g fixé et point stationnaire
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CORK = os.path.join(HERE, "..", "cork_bprime")

Y_CHAR = 0.20                      # rendement en char, les deux matériaux
K = (1.0 - Y_CHAR) / Y_CHAR        # B'g/B'c stationnaire = 4.0

MATERIALS = {
    "chêne": {
        "dir": HERE, "prefix": "oak",
        "rho_v": 700.0,            # ordre de grandeur (oak_bprime/README.md)
        "gas": (0.221, 0.536, 0.243),
        "color": "tab:brown",
    },
    "liège P50": {
        "dir": CORK, "prefix": "cork",
        "rho_v": 465.6,            # mesurée (P50)
        "gas": (0.287, 0.592, 0.121),
        "color": "tab:green",
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
    bc = read_csv(os.path.join(d, p + "_bprime_bc_table.csv"))
    ss = read_csv(os.path.join(d, p + "_bprime_steady_state.csv"))
    table = {}
    for r in bc:
        if abs(r["P_bar"] - P_BAR_1ATM) < 1e-6:
            table.setdefault(r["Bg"], {})[int(r["Tw_K"])] = r["Bc"]
    steady = {int(r["T_K"]): r for r in ss
              if abs(r["P_atm"] - 1.0) < 1e-9 and r["Bg_ss"] < BG_MAX_SWEEP}
    return table, steady


def main():
    data = {name: load(m) for name, m in MATERIALS.items()}

    print("Gaz de pyrolyse (fractions molaires élémentaires)")
    for name, m in MATERIALS.items():
        c, h, o = m["gas"]
        print(f"  {name:10s} C {c:.3f}  H {h:.3f}  O {o:.3f}  O/C = {o / c:.2f}")
    print(f"\nRendement en char {Y_CHAR:.0%} pour les deux : k = {K:.1f}\n")

    print("B'c à 1 atm, B'g imposé")
    print("  T [K] | " + " | ".join(
        f"B'g={bg:g} chêne / liège" for bg in BG_REPORT))
    for T in T_REPORT:
        cells = []
        for bg in BG_REPORT:
            vo = data["chêne"][0][bg][T]
            vc = data["liège P50"][0][bg][T]
            cells.append(f"{vo:.3f} / {vc:.3f}")
        print(f"  {T:5d} | " + " | ".join(cells))

    print("\nPoint stationnaire (B'g = 4 B'c), 1 atm")
    print("  T [K] | B'c chêne | B'c liège | B'g chêne | B'g liège | "
          "(1+k)B'c/ρv chêne | liège [1e-3 m3/kg] | rapport")
    rows = []
    for T in T_REPORT:
        so = data["chêne"][1].get(T)
        sc = data["liège P50"][1].get(T)
        if so is None or sc is None:
            continue
        ro = (1 + K) * so["Bc_ss"] / MATERIALS["chêne"]["rho_v"]
        rc = (1 + K) * sc["Bc_ss"] / MATERIALS["liège P50"]["rho_v"]
        print(f"  {T:5d} | {so['Bc_ss']:.4f}    | {sc['Bc_ss']:.4f}    | "
              f"{so['Bg_ss']:.3f}     | {sc['Bg_ss']:.3f}     | "
              f"{1e3 * ro:.3f}             | {1e3 * rc:.3f}              | "
              f"x{ro / rc:.2f}  (B'c x{so['Bc_ss'] / sc['Bc_ss']:.2f})")
        rows.append((T, so["Bc_Bg0"], so["Bc_ss"], sc["Bc_ss"],
                     so["Bg_ss"], sc["Bg_ss"], so["hw_ss_MJkg"], sc["hw_ss_MJkg"],
                     ro, rc))

    out_csv = os.path.join(HERE, "comparaison_chene_liege.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "Bc_Bg0", "Bc_ss_oak", "Bc_ss_corkP50",
                    "Bg_ss_oak", "Bg_ss_corkP50", "hw_ss_oak_MJkg",
                    "hw_ss_corkP50_MJkg", "recession_over_mdote_oak_m3kg",
                    "recession_over_mdote_corkP50_m3kg"])
        for r in rows:
            w.writerow([r[0]] + [f"{v:.6e}" for v in r[1:]])
    print(f"\n-> {os.path.relpath(out_csv, HERE)}")

    # --- figure ------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    styles = {0.0: ":", 0.5: "--", 2.0: "-.", 5.0: "-"}
    for name, m in MATERIALS.items():
        table, steady = data[name]
        for bg, ls in styles.items():
            if bg == 0.0 and name == "liège P50":
                continue            # identique au chêne à B'g = 0
            Ts = sorted(table[bg])
            lab = "B'g = 0 (les deux)" if bg == 0.0 else f"{name}, B'g = {bg:g}"
            ax1.plot(Ts, [table[bg][t] for t in Ts], ls,
                     color="k" if bg == 0.0 else m["color"], label=lab)
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
    fig.suptitle("Chêne vs liège/phénolique P50 — même char, même k, gaz opposés")
    fig.tight_layout()
    out_png = os.path.join(HERE, "comparaison_chene_liege.png")
    fig.savefig(out_png, dpi=150)
    print(f"-> {os.path.relpath(out_png, HERE)}")


if __name__ == "__main__":
    main()
