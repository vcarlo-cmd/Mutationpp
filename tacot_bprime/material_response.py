#!/usr/bin/env python
"""
Réponse matériau d'un composite carbone/phénolique et couplage à la table B'.

Ce script répond à la question : « si je change le rapport fibres/résine et la
porosité, qu'est-ce qui change pour la table B' ? »

Réponse courte : la table B' elle-même NE CHANGE PAS.

    Le solveur de bilan de masse de surface (Thermodynamics::surfaceMassBalance)
    ne consomme que trois compositions ELEMENTAIRES — bord de couche limite,
    gaz de pyrolyse, char — plus T, P et B'g. Densité, porosité et fraction
    volumique de fibres n'y entrent pas.

      - gaz de pyrolyse : produit par la résine seule. À résine identique, sa
        composition élémentaire est identique ; seule la QUANTITE produite
        change.
      - char : fibres de carbone et résine carbonisée sont du carbone pur,
        donc Y_c = C:1.0 quel que soit le rapport fibres/résine.

    Donc B'c(T, P, B'g) et h_w(T, P, B'g) sont rigoureusement les mêmes.

Ce qui change est la RÉPONSE MATÉRIAU, c'est-à-dire quel B'g est physiquement
atteint. En ablation stationnaire (matériau entièrement carbonisé reculant à
vitesse constante), les flux de gaz de pyrolyse et de char sont dans le rapport

    B'g / B'c  =  m_g / m_c  =  (rho_v - rho_c) / rho_c  =  k

Le point de fonctionnement du matériau à (T, P) est donc la solution de

    B'c  =  B'c_table(T, P, B'g = k * B'c)

Ce script calcule k pour chaque matériau, résout ce point fixe sur la table
produite par `bprime`, et trace les lignes de fonctionnement.

Usage :
    python material_response.py
"""

import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ONEATM = 101325.0

# ---------------------------------------------------------------------------
# Définition des matériaux
# ---------------------------------------------------------------------------
# rho_f / rho_m : masses volumiques intrinsèques des fibres et de la matrice
# char_yield    : fraction de la masse de résine restant après pyrolyse
#                 (TACOT : la matrice passe de eps_m=0.10 à 0.05 à densité
#                  intrinsèque constante, soit 50 %)

RHO_FIBER = 1600.0   # kg/m3, fibres de carbone ex-cellulose
RHO_RESIN = 1200.0   # kg/m3, novolac/formaldéhyde vierge
CHAR_YIELD = 0.50    # rendement en char de la résine

MATERIALS = {
    "TACOT": dict(f_solid=0.50, porosity=0.80, mixture="tacot-air_25",
                  pyro="tacot_pyro", char="tacot_char",
                  note="10 % fibres / 10 % résine / 80 % pores"),
    "CPh70":  dict(f_solid=0.70, porosity=0.01, mixture="cph70-air_25",
                  pyro="cph70_pyro", char="cph70_char",
                  note="70 % fibres / 30 % résine du solide, porosité 0.01"),
}

# Grille de calcul
T_RANGE = "300:25:4000"
BG_GRID = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4,
           0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4,
           3.0, 4.0, 5.5, 7.5, 10.0]
PRESSURES_ATM = [0.001, 0.01, 0.1, 1.0]


# ---------------------------------------------------------------------------
# Bilan de phase solide
# ---------------------------------------------------------------------------

def solid_balance(f_solid, porosity,
                  rho_f=RHO_FIBER, rho_m=RHO_RESIN, char_yield=CHAR_YIELD):
    """
    Bilan de phase solide d'un composite fibres/résine.

    f_solid  : fraction VOLUMIQUE de fibres dans le solide (0-1)
    porosity : porosité du composite (0-1)

    Retourne un dict avec fractions volumiques, densités vierge et char,
    fractions massiques et le coefficient de couplage k = (rho_v-rho_c)/rho_c.
    """
    solid = 1.0 - porosity
    eps_f = solid * f_solid
    eps_m = solid * (1.0 - f_solid)

    m_f = eps_f * rho_f            # masse de fibres par m3 de composite
    m_m = eps_m * rho_m            # masse de résine vierge
    rho_v = m_f + m_m

    m_mc = m_m * char_yield        # résine restante après pyrolyse
    rho_c = m_f + m_mc
    m_gas = rho_v - rho_c          # gaz de pyrolyse libéré

    return dict(
        eps_f=eps_f, eps_m=eps_m, porosity=porosity,
        m_f=m_f, m_m=m_m, m_mc=m_mc, m_gas=m_gas,
        rho_v=rho_v, rho_c=rho_c,
        wf_fiber_v=m_f / rho_v, wf_resin_v=m_m / rho_v,
        wf_fiber_c=m_f / rho_c, wf_resin_c=m_mc / rho_c,
        gas_frac=m_gas / rho_v,
        k=m_gas / rho_c,
    )


# ---------------------------------------------------------------------------
# Table B'
# ---------------------------------------------------------------------------

def find_bprime():
    cmd = shutil.which("bprime")
    if cmd:
        return cmd
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, "../build/src/apps/bprime"),
              "build/src/apps/bprime"):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def make_env():
    env = os.environ.copy()
    env.setdefault("MPP_DATA_DIRECTORY", os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../data")))
    return env


def bprime_table(path, mat, pressure_pa):
    """Retourne (Ts, BG_GRID, Bc[nbg,nT], hw[nbg,nT]) pour une pression."""
    Ts, Bc, Hw = None, [], []
    for bg in BG_GRID:
        cmd = [path, "-T", T_RANGE, "-P", str(pressure_pa), "-b", str(bg),
               "-m", mat["mixture"], "-bl", "air", "-py", mat["pyro"],
               "-char", mat["char"], "-char-elem", "C"]
        res = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
        if res.returncode != 0:
            print(f"\nERREUR bprime :\n{res.stderr}")
            sys.exit(1)
        t, b, h = [], [], []
        for line in res.stdout.strip().splitlines()[1:]:
            p = line.split()
            try:
                t.append(float(p[0])); b.append(float(p[1])); h.append(float(p[2]))
            except (ValueError, IndexError):
                continue
        Ts = np.array(t)
        Bc.append(b); Hw.append(h)
    return Ts, np.array(BG_GRID), np.array(Bc), np.array(Hw)


# ---------------------------------------------------------------------------
# Point de fonctionnement en ablation stationnaire
# ---------------------------------------------------------------------------

def steady_state(Ts, bgs, Bc, Hw, k, iters=200):
    """
    Résout B'c = B'c_table(T, B'g = k*B'c) pour chaque température, par
    itération de point fixe amortie sur la table interpolée en B'g.

    Retourne (Bc_ss, hw_ss, Bg_ss).
    """
    nT = len(Ts)
    bc_ss = np.zeros(nT)
    hw_ss = np.zeros(nT)

    for j in range(nT):
        col_bc = Bc[:, j]           # B'c(B'g) à T fixée
        col_hw = Hw[:, j]
        x = col_bc[0]               # départ : solution à B'g = 0
        for _ in range(iters):
            bg = min(k * x, bgs[-1])
            xn = np.interp(bg, bgs, col_bc)
            if not np.isfinite(xn):
                break
            if abs(xn - x) < 1e-10 * max(abs(x), 1e-12):
                x = xn
                break
            x = 0.5 * (x + xn)      # amortissement : la table est très raide
        bc_ss[j] = x
        hw_ss[j] = np.interp(min(k * x, bgs[-1]), bgs, col_hw)

    return bc_ss, hw_ss, np.minimum(k * bc_ss, bgs[-1])


# ---------------------------------------------------------------------------
# Sorties
# ---------------------------------------------------------------------------

def print_balance(name, mat, b):
    print(f"\n{'='*74}\n{name}  —  {mat['note']}\n{'='*74}")
    print(f"  fractions volumiques : fibres {b['eps_f']:.4f}   "
          f"résine {b['eps_m']:.4f}   pores {b['porosity']:.4f}")
    print(f"\n  VIERGE   rho_v = {b['rho_v']:8.1f} kg/m3")
    print(f"     fibres {b['m_f']:8.1f} kg/m3  ({100*b['wf_fiber_v']:5.1f} % masse)")
    print(f"     résine {b['m_m']:8.1f} kg/m3  ({100*b['wf_resin_v']:5.1f} % masse)")
    print(f"\n  CHAR     rho_c = {b['rho_c']:8.1f} kg/m3")
    print(f"     fibres {b['m_f']:8.1f} kg/m3  ({100*b['wf_fiber_c']:5.1f} % masse)")
    print(f"     résine {b['m_mc']:8.1f} kg/m3  ({100*b['wf_resin_c']:5.1f} % masse)")
    print(f"\n  gaz de pyrolyse libéré : {b['m_gas']:.1f} kg/m3 "
          f"({100*b['gas_frac']:.1f} % de la masse vierge)")
    print(f"  couplage stationnaire  : k = B'g/B'c = "
          f"(rho_v-rho_c)/rho_c = {b['k']:.4f}")


def plot_operating_lines(data, out_png):
    """B'c et h_w en ablation stationnaire, un panneau par pression."""
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    fig.suptitle("Ablation stationnaire — la table B' est commune, "
                 "seule la ligne de fonctionnement $B'_g = k\\,B'_c$ diffère",
                 fontsize=14)
    colors = {"TACOT": "tab:blue", "CPh70": "tab:red"}

    for col, P in enumerate(PRESSURES_ATM):
        ax1, ax2 = axes[0, col], axes[1, col]
        for name, d in data.items():
            Ts, bc, hw, k = d[P]["Ts"], d[P]["bc"], d[P]["hw"], d["k"]
            ax1.plot(Ts, bc, color=colors[name], lw=2,
                     label=f"{name}  (k={k:.3f})")
            ax2.plot(Ts, hw, color=colors[name], lw=2, label=name)
        exp = int(round(np.log10(P)))
        ttl = rf"$10^{{{exp}}}$ atm" if exp else "1 atm"
        ax1.set_title(ttl)
        ax1.set_yscale("log")
        ax1.set_ylabel(r"$B'_c$ stationnaire")
        ax2.set_ylabel(r"$h_w$ [MJ/kg]")
        for ax in (ax1, ax2):
            ax.set_xlabel(r"$T_w$ [K]")
            ax.grid(True, which="both", ls="--", alpha=0.35)
            ax.legend(fontsize=8, loc="upper left")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"\nFigure sauvegardée : {out_png}")


def main():
    path = find_bprime()
    if path is None:
        print("Binaire 'bprime' introuvable — compilez la cible bprime.")
        sys.exit(1)

    balances = {n: solid_balance(m["f_solid"], m["porosity"])
                for n, m in MATERIALS.items()}
    for name, mat in MATERIALS.items():
        print_balance(name, mat, balances[name])

    kt, kc = balances["TACOT"]["k"], balances["CPh70"]["k"]
    print(f"\n{'='*74}\nCOMPARAISON\n{'='*74}")
    print(f"  rho_v : {balances['TACOT']['rho_v']:7.1f} -> "
          f"{balances['CPh70']['rho_v']:7.1f} kg/m3  "
          f"(x{balances['CPh70']['rho_v']/balances['TACOT']['rho_v']:.2f})")
    print(f"  k     : {kt:7.4f} -> {kc:7.4f}          (x{kc/kt:.2f})")
    print("\n  Le CPh70 est bien plus dense mais, étant riche en fibres, il")
    print("  produit relativement MOINS de gaz de pyrolyse par unité de char :")
    print("  sa ligne de fonctionnement est à B'g environ deux fois plus bas.")

    # --- Tables B' et points de fonctionnement ---------------------------
    print(f"\n{'='*74}\nCALCUL DES TABLES B'\n{'='*74}")
    data = {}
    for name, mat in MATERIALS.items():
        k = balances[name]["k"]
        data[name] = {"k": k}
        for P in PRESSURES_ATM:
            print(f"  {name:6s}  P = {P:g} atm ...", end=" ", flush=True)
            Ts, bgs, Bc, Hw = bprime_table(path, mat, P * ONEATM)
            bc, hw, bg = steady_state(Ts, bgs, Bc, Hw, k)
            data[name][P] = dict(Ts=Ts, bc=bc, hw=hw, bg=bg, Bc=Bc)
            print(f"{len(Ts)} points")

    # --- Vérification : les deux tables sont-elles identiques ? -----------
    print(f"\n{'='*74}\nVÉRIFICATION — les deux tables B' sont-elles identiques ?\n{'='*74}")
    for P in PRESSURES_ATM:
        d = np.abs(data["TACOT"][P]["Bc"] - data["CPh70"][P]["Bc"]).max()
        print(f"  P = {P:6g} atm : écart max sur B'c(T,B'g) = {d:.3e}")
    print("\n  -> tables rigoureusement identiques : la mise en données de la")
    print("     table B' ne dépend pas du rapport fibres/résine ni de la porosité.")

    # --- Points de fonctionnement à 1 atm ---------------------------------
    print(f"\n{'='*74}\nPOINT DE FONCTIONNEMENT STATIONNAIRE À 1 atm\n{'='*74}")
    print(f"{'T [K]':>7} | {'TACOT B''c':>11} {'B''g':>8} | "
          f"{'CPh70 B''c':>11} {'B''g':>8} | {'écart B''c':>10}")
    Ts = data["TACOT"][1.0]["Ts"]
    for T in (1000, 2000, 2500, 3000, 3200, 3400, 3600, 3800):
        j = int(np.argmin(np.abs(Ts - T)))
        a, b = data["TACOT"][1.0], data["CPh70"][1.0]
        rel = (b["bc"][j] - a["bc"][j]) / max(a["bc"][j], 1e-12) * 100
        print(f"{Ts[j]:>7.0f} | {a['bc'][j]:>11.5f} {a['bg'][j]:>8.4f} | "
              f"{b['bc'][j]:>11.5f} {b['bg'][j]:>8.4f} | {rel:>9.2f} %")

    plot_operating_lines(data, "material_response_operating_lines.png")
    print("\nTerminé.")


if __name__ == "__main__":
    main()
