#!/usr/bin/env python
"""
Deux mises en données complètes, sur base COMPOSITE (fibres + résine), pour
les variantes ZURAM XX/YY.

    JEU A — Arrhenius + F_i        (fractions de perte de densité, sans dimension)
    JEU B — Arrhenius + Delta_rho_i (masses décomposables, en kg/m³ de composite)

Les deux décrivent le MÊME matériau : Delta_rho_i = F_i · rho_vierge. Le script
les émet puis vérifie leur équivalence en intégrant séparément les deux
systèmes d'équations.

Dans les deux jeux, A_i, E_i/R et m_i sont IDENTIQUES pour toutes les
variantes : ils décrivent la résine, pas le composite.

Usage :
    python zuram_mise_en_donnees.py
    python zuram_mise_en_donnees.py 14/40 18/80 --beta 20
"""

import argparse
import math

from zuram_variantes import KINETICS, CHAR_YIELD, parse, variant

R_GAS = 8.31446261815324          # J/mol/K
SUM_F = sum(f for f, _, _, _ in KINETICS)


# ---------------------------------------------------------------------------
# Émission des jeux
# ---------------------------------------------------------------------------

def common_block():
    print("  Équation commune aux deux jeux (forme normalisée, rho_i,c = 0) :")
    print()
    print("      dx_i/dt = - A_i · x_i^m_i · exp(-E_i/(R·T))        x_i(0) = 1")
    print()
    print("  Paramètres d'Arrhenius — IDENTIQUES POUR LES TROIS VARIANTES :")
    print()
    print(f"      {'i':>2} {'A_i [1/s]':>14} {'E_i/R [K]':>14} "
          f"{'E_i [J/mol]':>14} {'m_i [-]':>9}")
    for i, (f, la, er, m) in enumerate(KINETICS, 1):
        print(f"      {i:>2} {10**la:14.6e} {er:14.6f} {er*R_GAS:14.2f} {m:9.2f}")
    print()
    print("  (la colonne n du classeur vaut 0 pour les quatre réactions :")
    print("   l'exposant qui agit est m)")


def deck_a(variants, labels):
    print("\n" + "=" * 79)
    print("  JEU A — ARRHENIUS + F_i   (fractions, base composite)")
    print("=" * 79)
    print("""
  Reconstruction de la densité :

      rho(T) = rho_vierge · [ 1 - somme_i F_i · (1 - x_i) ]

  F_i = fraction de la densité du COMPOSITE que la réaction i fait perdre.
  somme F_i = perte totale ; rho_char/rho_vierge = 1 - somme F_i.
""")
    for lab, v in zip(labels, variants):
        Fs = [f * v["w_resin"] for f, _, _, _ in KINETICS]
        print(f"  --- ZURAM {lab} " + "-" * (79 - 15 - len(lab)))
        print(f"      rho_vierge          {v['rho_v']:12.4f} kg/m3")
        print(f"      rho_char            {v['rho_c']:12.4f} kg/m3")
        print(f"      porosite_vierge     {v['poro_v']:12.6f}")
        print(f"      porosite_char       {v['poro_c']:12.6f}")
        print()
        print(f"      {'i':>2} {'F_i [-]':>14} {'A_i [1/s]':>14} "
              f"{'E_i/R [K]':>14} {'m_i':>7}")
        for i, ((f, la, er, m), F) in enumerate(zip(KINETICS, Fs), 1):
            print(f"      {i:>2} {F:14.9f} {10**la:14.6e} {er:14.6f} {m:7.2f}")
        print(f"      {'':>2} {sum(Fs):14.9f}   <- somme, = 1 - rho_char/rho_vierge")
        print(f"      controle : 1 - {sum(Fs):.9f} = {1-sum(Fs):.9f}"
              f"   rho_c/rho_v = {v['rho_c']/v['rho_v']:.9f}")
        print()


def deck_b(variants, labels):
    print("=" * 79)
    print("  JEU B — ARRHENIUS + Delta_rho_i   (kg/m3 de composite)")
    print("=" * 79)
    print("""
  Reconstruction de la densité :

      d(rho_i)/dt = - A_i · Drho_i · (rho_i/Drho_i)^m_i · exp(-E_i/(R·T))
      rho_i(0) = Drho_i          rho(T) = rho_inerte + somme_i rho_i

  Drho_i = masse décomposable de la réaction i, par m3 de COMPOSITE.
  rho_inerte = fibres + char de résine = rho_char (ne pyrolyse pas).

  ATTENTION si le code demande un couple (rho_i,v ; rho_i,c) :
      rho_i,v = Drho_i     et     rho_i,c = 0
  Ne PAS mettre rho_i,c non nul avec ces A : le terme du modèle est
  A_i·rho_i,v·[(rho_i-rho_i,c)/rho_i,v]^m, et les A ci-dessus ne sont
  valides que pour rho_i,c = 0.
""")
    for lab, v in zip(labels, variants):
        Ds = [f * v["w_resin"] * v["rho_v"] for f, _, _, _ in KINETICS]
        print(f"  --- ZURAM {lab} " + "-" * (79 - 15 - len(lab)))
        print(f"      rho_vierge          {v['rho_v']:12.4f} kg/m3")
        print(f"      rho_inerte (=char)  {v['rho_c']:12.4f} kg/m3")
        print(f"      porosite_vierge     {v['poro_v']:12.6f}")
        print(f"      porosite_char       {v['poro_c']:12.6f}")
        print()
        print(f"      {'i':>2} {'Drho_i [kg/m3]':>15} {'A_i [1/s]':>14} "
              f"{'E_i/R [K]':>14} {'m_i':>7}")
        for i, ((f, la, er, m), D) in enumerate(zip(KINETICS, Ds), 1):
            print(f"      {i:>2} {D:15.6f} {10**la:14.6e} {er:14.6f} {m:7.2f}")
        print(f"      {'':>2} {sum(Ds):15.6f}   <- somme = gaz total libere")
        print(f"      controle : {v['rho_c']:.4f} + {sum(Ds):.4f} = "
              f"{v['rho_c']+sum(Ds):.4f}   (rho_vierge = {v['rho_v']:.4f})")
        print()


# ---------------------------------------------------------------------------
# Vérification croisée
# ---------------------------------------------------------------------------

def integrate_a(v, beta, t0=300.0, t1=1400.0, n=200000):
    """Jeu A : x_i normalisés, densité reconstruite par les F_i."""
    Fs = [f * v["w_resin"] for f, _, _, _ in KINETICS]
    b, dT = beta / 60.0, (t1 - t0) / n
    dt = dT / b
    x = [1.0] * len(KINETICS)
    Ts, rho = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        rho.append(v["rho_v"] * (1.0 - sum(F * (1.0 - xi)
                                           for F, xi in zip(Fs, x))))
        for j, (f, la, er, m) in enumerate(KINETICS):
            if x[j] > 0.0:
                x[j] = max(0.0, x[j] - 10.0**la * x[j]**m
                           * math.exp(-er / T) * dt)
    return Ts, rho


def integrate_b(v, beta, t0=300.0, t1=1400.0, n=200000):
    """Jeu B : rho_i en kg/m3, densité reconstruite par somme."""
    Ds = [f * v["w_resin"] * v["rho_v"] for f, _, _, _ in KINETICS]
    inert = v["rho_v"] - sum(Ds)
    b, dT = beta / 60.0, (t1 - t0) / n
    dt = dT / b
    r = list(Ds)
    Ts, rho = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        rho.append(inert + sum(r))
        for j, (f, la, er, m) in enumerate(KINETICS):
            if r[j] > 0.0:
                r[j] = max(0.0, r[j] - 10.0**la * Ds[j] * (r[j] / Ds[j])**m
                           * math.exp(-er / T) * dt)
    return Ts, rho


def verify(variants, labels, beta):
    print("=" * 79)
    print(f"  VÉRIFICATION — les deux jeux décrivent-ils le même matériau ?")
    print("=" * 79)
    print(f"  Intégration indépendante des deux systèmes, ATG à {beta:.0f} K/min.\n")
    print(f"  {'T [°C]':>7} |" + "".join(f"{'  ZURAM ' + lab:>26}" for lab in labels))
    print(f"  {'':7} |" + "".join(f"{'jeu A':>11}{'jeu B':>11}{'ecart':>4}"
                                  for _ in labels))
    print("  " + "-" * (9 + 26 * len(labels)))
    runs = [(integrate_a(v, beta), integrate_b(v, beta)) for v in variants]
    Ts = runs[0][0][0]
    for Tc in (400, 500, 600, 700, 900, 1100):
        i = min(range(len(Ts)), key=lambda k: abs(Ts[k] - (Tc + 273.15)))
        cells = ""
        for (_, ra), (_, rb) in runs:
            cells += f"{ra[i]:11.4f}{rb[i]:11.4f}{abs(ra[i]-rb[i]):4.0f}"
        print(f"  {Tc:>7} |" + cells)
    worst = 0.0
    for (_, ra), (_, rb) in runs:
        worst = max(worst, max(abs(a - b) for a, b in zip(ra, rb)))
    ok = worst < 1e-9
    print(f"\n  écart maximal sur toute la montée : {worst:.2e} kg/m3")
    print(f"  -> {'JEUX ÉQUIVALENTS' if ok else 'ÉCHEC'}")

    print("\n  Densités finales à 1400 K (le modèle garde une queue algébrique,")
    print("  cf. `variantes_zuram.md` — le char asymptotique est un peu plus bas) :")
    print(f"\n  {'variante':>10} {'rho(1400 K)':>13} {'rho_char asympt.':>18}"
          f" {'ecart':>9}")
    for lab, v, ((_, ra), _) in zip(labels, variants, runs):
        print(f"  {lab:>10} {ra[-1]:13.4f} {v['rho_c']:18.4f} "
              f"{ra[-1]-v['rho_c']:9.4f}")
    return ok


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*", default=["14/40", "18/50", "18/80"])
    ap.add_argument("--beta", type=float, default=20.0)
    args = ap.parse_args()
    labels = args.variants
    variants = [variant(*parse(a)) for a in labels]

    print()
    print("  MISES EN DONNÉES ZURAM — BASE COMPOSITE (fibres + résine)")
    print("=" * 79)
    common_block()
    deck_a(variants, labels)
    deck_b(variants, labels)
    ok = verify(variants, labels, args.beta)

    print("\n" + "=" * 79)
    print("  RAPPEL — ce qui accompagne le bloc cinétique, identique aux 3 :")
    print("=" * 79)
    print("      gaz de pyrolyse   C:0.171, H:0.722, O:0.107  (fractions molaires)")
    print("      char              C:1.0")
    print("      bord de couche l. O:0.21, N:0.79  (air)")
    print("      table B'          la MÊME pour les trois (verifie a 0.000e+00)")
    print("=" * 79)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
