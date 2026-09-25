#!/usr/bin/env python
"""
Table B' du liège/phénolique (cork P50) sous jet de combustion kérosène
Jet-A1 / oxygène — torche HVOF (type JP-5000) — et comparaison à la table
air. Cas distinct de l'OAT (cork_bprime_oat.py).

Seul le bord de couche limite change par rapport à cork_bprime.py : même
gaz de pyrolyse (cork_pyro), même char (C pur), même liste d'espèces
(data/mixtures/cork-hvof.xml reprend celle de cork-air.xml).

Carburant : Jet-A1 = C12H23 (formule moyenne, M = 167.3 g/mol).
Comburant : O2 pur. Richesse phi :

    C12H23 + (17.75/phi) O2          éléments C:H:O = 12 : 23 : 35.5/phi

    phi = 0.8  pauvre     phi = 1.0  stoechiométrique     phi = 1.3  riche

Jet + air entraîné : une fraction MASSIQUE f d'air ambiant se mélange aux
produits avant la paroi (f = 0.25 et 0.50, pour chaque phi). L'azote entre
alors au bord.

Plateau d'oxydation attendu (tout l'O en excès sur C forme du CO avec le
carbone du char, H finit en H2) :

    B'c = M_C (35.5/phi - 12) / (M_C12H23 + (17.75/phi) M_O2)
        phi = 0.8 : 0.443    phi = 1.0 : 0.384    phi = 1.3 : 0.304
        (air      : 0.175)

et, avec entraînement, la moyenne massique (1 - f) B'c(jet) + f B'c(air).
Contrairement à l'OAT neutre, le jet HVOF est bien plus oxydant que l'air :
l'air entraîné DIMINUE ici B'c.

Sorties :
    cork_hvof_bprime_bc_table.csv     Tw_K, P_bar, phi, f_air, Bg, Bc
    cork_hvof_bprime_hw_table.csv     Tw_K, P_bar, phi, f_air, Bg, hw_Jkg
    cork_hvof_bprime_steady_state.csv point de fonctionnement k = 4 (+ air)
    cork_hvof_bprime_vs_air.png       B'c et h_w à 1 atm, B'g = 0 et 5
    cork_hvof_bprime_steady_state.png B'c stationnaire, air vs HVOF
    cork_hvof_air_bprime.png          effet de l'air entraîné (1 atm)

Usage :
    python cork_bprime_hvof.py
"""

import csv
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cork_bprime import (find_bprime, make_env, parse_output,
                         solve_operating_point, K_CORK, RHO_CHAR, ONEATM)

MIXTURE   = "cork-hvof"
PYRO_COMP = "cork_pyro"
CHAR_COMP = "cork_char"
CHAR_ELEM = "C"

# (étiquette composition XML, richesse phi, fraction massique d'air, libellé)
# phi = None : bord de couche limite air (référence).
EDGES = [
    ("air",               None, 1.00, "air"),
    ("hvof_phi0p8",       0.8,  0.00, r"HVOF $\phi$ = 0.8 (pauvre)"),
    ("hvof_phi1p0",       1.0,  0.00, r"HVOF $\phi$ = 1.0 (stoech.)"),
    ("hvof_phi1p3",       1.3,  0.00, r"HVOF $\phi$ = 1.3 (riche)"),
    ("hvof_phi0p8_air25", 0.8,  0.25, r"$\phi$ = 0.8 + 25 % air"),
    ("hvof_phi0p8_air50", 0.8,  0.50, r"$\phi$ = 0.8 + 50 % air"),
    ("hvof_phi1p0_air25", 1.0,  0.25, r"$\phi$ = 1.0 + 25 % air"),
    ("hvof_phi1p0_air50", 1.0,  0.50, r"$\phi$ = 1.0 + 50 % air"),
    ("hvof_phi1p3_air25", 1.3,  0.25, r"$\phi$ = 1.3 + 25 % air"),
    ("hvof_phi1p3_air50", 1.3,  0.50, r"$\phi$ = 1.3 + 50 % air"),
]
PHI   = {bl: phi for bl, phi, _, _ in EDGES}
F_AIR = {bl: f for bl, _, f, _ in EDGES}
LABEL = {bl: lab for bl, _, _, lab in EDGES}
PURE  = [bl for bl, phi, f, _ in EDGES if phi is None or f == 0.0]

T_RANGE       = "300:25:5000"
PRESSURES_ATM = np.logspace(-3, 3, 25)
BG_VALUES     = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

BG_SWEEP = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4,
            0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4,
            3.0, 4.0, 5.5, 7.5, 10.0]
SWEEP_T_RANGE   = "300:25:4000"
SWEEP_PRESSURES = [0.01, 1.0]

M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007
N_C, N_H = 12, 23                       # Jet-A1 = C12H23
O2_ST = N_C + N_H / 4.0                 # 17.75
B_AIR = M_C * 0.21 / (0.79 * M_N + 0.21 * M_O)   # plateau de l'air, 0.175

COLORS = {0.8: "tab:green", 1.0: "tab:red", 1.3: "tab:purple", None: "k"}
STYLE  = {0.0: "-", 0.25: "--", 0.50: ":", 1.0: "-"}
T_REPORT = (1000, 1500, 2000, 2500, 3000, 3400)


def plateau(phi, f_air=0.0):
    """B'c de plateau d'oxydation (O en excès -> CO) pour C12H23 + O2 à la
    richesse phi, mélangé à une fraction massique f_air d'air entraîné."""
    n_o2 = O2_ST / phi
    m = N_C * M_C + N_H * M_H + 2 * n_o2 * M_O
    b = M_C * max(2 * n_o2 - N_C, 0) / m
    return (1.0 - f_air) * b + f_air * B_AIR


def run(bprime_path, bl, P_pa, bg, t_range=T_RANGE):
    cmd = [bprime_path, "-T", t_range, "-P", str(P_pa), "-b", str(bg),
           "-m", MIXTURE, "-bl", bl, "-py", PYRO_COMP,
           "-char", CHAR_COMP, "-char-elem", CHAR_ELEM]
    res = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if res.returncode != 0:
        print(f"ERREUR bprime ({bl}, P={P_pa} Pa, B'g={bg}) :\n{res.stderr}")
        sys.exit(1)
    return parse_output(res.stdout)[1]


def main():
    bprime_path = find_bprime()
    if bprime_path is None:
        sys.exit("Binaire 'bprime' introuvable (cf. cork_bprime.py).")

    # --- Tables B'c / h_w, format long ------------------------------------
    tables = {}   # (bl, P_atm, bg) -> data
    for bl, _, _, _ in EDGES:
        for bg in BG_VALUES:
            print(f"  {bl:17s} B'g = {bg:<4g}", flush=True)
            for P in PRESSURES_ATM:
                tables[(bl, P, bg)] = run(bprime_path, bl, P * ONEATM, bg)

    with open("cork_hvof_bprime_bc_table.csv", "w", newline="") as fb, \
         open("cork_hvof_bprime_hw_table.csv", "w", newline="") as fh:
        wb, wh = csv.writer(fb), csv.writer(fh)
        wb.writerow(["Tw_K", "P_bar", "phi", "f_air", "Bg", "Bc"])
        wh.writerow(["Tw_K", "P_bar", "phi", "f_air", "Bg", "hw_Jkg"])
        for bl, phi, f, _ in EDGES:
            if phi is None:
                continue   # la table air est déjà dans cork_bprime_*.csv
            for bg in BG_VALUES:
                for P in PRESSURES_ATM:
                    P_bar = P * ONEATM / 1.0e5
                    for row in tables[(bl, P, bg)]:
                        key = [f"{row[0]:.6g}", f"{P_bar:.6g}", f"{phi:g}",
                               f"{f:g}", f"{bg:g}"]
                        wb.writerow(key + [f"{row[1]:.6e}"])
                        wh.writerow(key + [f"{row[2] * 1e6:.6e}"])
    print("Tables : cork_hvof_bprime_bc_table.csv, "
          "cork_hvof_bprime_hw_table.csv")

    # --- Figure B'c / h_w à 1 atm, jet pur vs air -------------------------
    P1 = PRESSURES_ATM[np.argmin(np.abs(np.log10(PRESSURES_ATM)))]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Table B' — liège/phénolique P50 : air vs jet HVOF "
                 "Jet-A1/O$_2$ (1 atm)", fontsize=13)
    for bl in PURE:
        c = COLORS[PHI[bl]]
        for bg, ls in ((0.0, "-"), (5.0, "--")):
            d = tables[(bl, P1, bg)]
            ax1.plot(d[:, 0], np.maximum(d[:, 1], 1e-5), ls, color=c, lw=2,
                     label=f"{LABEL[bl]}, B'g = {bg:g}")
            ax2.plot(d[:, 0], d[:, 2], ls, color=c, lw=2)
        if PHI[bl] is not None:
            ax1.axhline(plateau(PHI[bl]), color=c, lw=0.8, ls=":")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-4, 3e2)
    ax1.set_xlabel(r"$T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$")
    ax1.set_title(r"$B'_c$ (pointillés : plateau O en excès $\to$ CO)")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=7.5, loc="upper left")
    ax2.set_xlabel(r"$T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"$h_w$ (trait plein B'g = 0, tirets B'g = 5)")
    ax2.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig("cork_hvof_bprime_vs_air.png", dpi=150)
    plt.close()

    # --- Point de fonctionnement stationnaire (k = 4) ---------------------
    ss = {}
    for bl, _, _, _ in EDGES:
        for P in SWEEP_PRESSURES:
            raw = {bg: run(bprime_path, bl, P * ONEATM, bg, SWEEP_T_RANGE)
                   for bg in BG_SWEEP}
            Ts = raw[0.0][:, 0]
            rows = []
            for j, T in enumerate(Ts):
                bcs = np.array([raw[bg][j, 1] for bg in BG_SWEEP])
                hws = np.array([raw[bg][j, 2] for bg in BG_SWEEP])
                bc = solve_operating_point(np.array(BG_SWEEP), bcs, K_CORK)
                bg = min(K_CORK * bc, BG_SWEEP[-1])
                rows.append((T, bc, bg, np.interp(bg, BG_SWEEP, hws), bcs[0]))
            ss[(bl, P)] = np.array(rows)
            print(f"  point de fonctionnement {bl:17s} P = {P:g} atm")

    with open("cork_hvof_bprime_steady_state.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["edge", "phi", "f_air", "P_atm", "T_K", "Bc_ss",
                    "Bg_ss", "hw_ss_MJkg", "Bc_Bg0",
                    "recession_over_mdote_m3_per_kg"])
        for bl, phi, fa, _ in EDGES:
            for P in SWEEP_PRESSURES:
                for T, bc, bg, hw, bc0 in ss[(bl, P)]:
                    w.writerow([bl, "" if phi is None else f"{phi:g}",
                                f"{fa:g}", f"{P:g}", f"{T:g}", f"{bc:.6e}",
                                f"{bg:.6e}", f"{hw:.6e}", f"{bc0:.6e}",
                                f"{bc / RHO_CHAR:.6e}"])
    print("Point de fonctionnement : cork_hvof_bprime_steady_state.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(r"Liège/phénolique P50 — $B'_c$ stationnaire "
                 r"($B'_g = 4\,B'_c$) : air vs jet HVOF", fontsize=13)
    for ax, P in zip(axes, SWEEP_PRESSURES):
        for bl in PURE:
            d = ss[(bl, P)]
            ok = d[:, 2] < BG_SWEEP[-1]
            ax.plot(d[ok, 0], np.maximum(d[ok, 1], 1e-5),
                    color=COLORS[PHI[bl]], lw=2, label=LABEL[bl])
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 3)
        ax.set_xlabel(r"$T_w$ [K]")
        ax.set_ylabel(r"$B'_c$ stationnaire")
        ax.set_title("1 atm" if P == 1.0 else f"{P:g} atm")
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=9, loc="upper left")
    plt.tight_layout()
    plt.savefig("cork_hvof_bprime_steady_state.png", dpi=150)
    plt.close()

    # --- Effet de l'air entraîné (1 atm) ----------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle("Liège/phénolique P50 — jet HVOF Jet-A1/O$_2$ + air "
                 "entraîné (fraction massique f), 1 atm", fontsize=13)
    for bl, phi, fa, lab in EDGES:
        c, ls = COLORS[phi], STYLE[fa]
        if phi is not None and fa == 0.0:
            lab = rf"$\phi$ = {phi:g} (sans air)"
        d = tables[(bl, P1, 0.0)]
        axes[0].plot(d[:, 0], np.maximum(d[:, 1], 1e-5), ls, color=c, lw=2,
                     label=lab)
        axes[1].plot(d[:, 0], d[:, 2], ls, color=c, lw=2, label=lab)
        d = ss[(bl, 1.0)]
        ok = d[:, 2] < BG_SWEEP[-1]
        axes[2].plot(d[ok, 0], np.maximum(d[ok, 1], 1e-5), ls, color=c,
                     lw=2, label=lab)
    for ax, yl, t in ((axes[0], r"$B'_c$", r"$B'_c$, $B'_g = 0$"),
                      (axes[1], r"$h_w$ [MJ/kg]", r"$h_w$, $B'_g = 0$"),
                      (axes[2], r"$B'_c$ stationnaire",
                       r"$B'_c$ stationnaire ($B'_g = 4\,B'_c$)")):
        ax.set_xlabel(r"$T_w$ [K]")
        ax.set_ylabel(yl)
        ax.set_title(t)
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=7.5, loc="upper left")
    for ax in (axes[0], axes[2]):
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 3)
    axes[0].set_xlim(300, 4000)
    axes[1].set_xlim(300, 4000)
    plt.tight_layout()
    plt.savefig("cork_hvof_air_bprime.png", dpi=150)
    plt.close()

    # --- Récapitulatif ----------------------------------------------------
    print("\nPlateaux théoriques : " + ", ".join(
        f"{bl} -> {plateau(PHI[bl], F_AIR[bl]):.4f}"
        for bl, _, _, _ in EDGES if PHI[bl]))
    print(f"\n1 atm   {'':17s}" + "".join(f"{T:>10d}" for T in T_REPORT))
    for kind, bg in (("B'c B'g=0", 0.0), ("B'c B'g=5", 5.0)):
        for bl, _, _, _ in EDGES:
            d = tables[(bl, P1, bg)]
            vals = [np.interp(T, d[:, 0], d[:, 1]) for T in T_REPORT]
            print(f"{kind} {bl:17s}" + "".join(f"{v:>10.4f}" for v in vals))
    for bl, _, _, _ in EDGES:
        d = ss[(bl, 1.0)]
        vals = [np.interp(T, d[:, 0], d[:, 1]) for T in T_REPORT]
        hws = [np.interp(T, d[:, 0], d[:, 3]) for T in T_REPORT]
        print(f"B'c ss    {bl:17s}" + "".join(f"{v:>10.4f}" for v in vals))
        print(f"hw ss     {bl:17s}" + "".join(f"{v:>10.2f}" for v in hws))


if __name__ == "__main__":
    main()
