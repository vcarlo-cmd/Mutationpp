#!/usr/bin/env python
"""
Mises en données de pyrolyse, dans les trois formes usuelles.

Traite les matériaux qui partagent les fibres et la résine du TACOT — le TACOT
lui-même et le CPh70 — puisque seule leur microstructure diffère : même
cinétique, seul eps_m change.

Pendant de `../zuram_bprime/zuram_mise_en_donnees.py`. Le TACOT pose une
difficulté que le ZURAM n'avait pas : sa phase B a une densité résiduelle
rho_c NON NULLE (900 -> 600). Les formes « en f » et « en Delta_rho », qui
supposent implicitement rho_i,c = 0, ne peuvent donc PAS reprendre les A du
classeur tels quels — il faut les corriger.

    JEU 0 — natif : (rho_i,v ; rho_i,c) + A du classeur    <- ce que publie TACOT_3.0
    JEU A — Arrhenius + F_i          (fractions, base composite)
    JEU B — Arrhenius + Delta_rho_i  (kg/m³ de composite)

Les jeux A et B exigent A_B / 9. Le script établit ce facteur, émet les trois
jeux, puis vérifie leur équivalence en intégrant les trois systèmes
séparément — et montre ce que coûte l'oubli de la correction.

Source : TACOT_3.0.xls, feuille 'Pyrolysis model' (Goldstein 1965/1969),
cellules A11:G20 et A31:A34. Voir `resine_tacot.md`.

Usage :
    python tacot_mise_en_donnees.py            # TACOT (défaut)
    python tacot_mise_en_donnees.py CPh70
    python tacot_mise_en_donnees.py CPh70 --beta 10
"""

import argparse
import math

# --- 'Pyrolysis model'!B18:G19, base RÉSINE (kg/m³ de matrice) --------------
#     rho_i_v, rho_i_c, A [1/s], E/R [K], psi, T_reac [K]
GOLDSTEIN = [
    ("A", 300.0,   0.0, 1.2e4,  8555.555555555555,  3.0, 333.3333333333333),
    ("B", 900.0, 600.0, 4.48e9, 20444.444444444445, 3.0, 555.5555555555555),
]
RHO_RESIN = 1200.0        # somme des rho_i_v = densité intrinsèque de la matrice
RHO_FIBER = 1600.0        # 'Pyrolysis model'!B20 — les fibres ne pyrolysent pas

PYRO_GAS = {"C": 0.206, "H": 0.679, "O": 0.115}   # 'Pyrolysis model'!B4:B6

# Matériaux partageant CES fibres et CETTE résine : seules les proportions et
# la porosité changent, donc seul eps_m change. La cinétique est la même.
#   f_solid  : fraction VOLUMIQUE de fibres dans le solide
#   porosity : porosité du composite vierge
MATERIALS = {
    "TACOT": dict(f_solid=0.50, porosity=0.80,
                  note="10 % fibres / 10 % résine / 80 % pores "
                       "('Description'!A15:A17)"),
    "CPh70": dict(f_solid=0.70, porosity=0.01,
                  note="70 % fibres / 30 % résine du solide, porosité 0.01 "
                       "(cf. ../tacot_bprime/autre_materiau.md)"),
}


def microstructure(name):
    """Fractions volumiques et densités du composite, depuis f_solid/porosity."""
    m = MATERIALS[name]
    solid = 1.0 - m["porosity"]
    eps_f = m["f_solid"] * solid
    eps_m = (1.0 - m["f_solid"]) * solid
    rho_v = eps_f * RHO_FIBER + eps_m * RHO_RESIN
    # char : la matrice perd somme(rho_i_v - rho_i_c) par m3 de matrice
    dmat = sum(g[1] - g[2] for g in GOLDSTEIN)
    rho_c = rho_v - eps_m * dmat
    return dict(name=name, note=m["note"], eps_f=eps_f, eps_m=eps_m,
                rho_v=rho_v, rho_c=rho_c, porosity=m["porosity"])


def corrected_A(rho_v, rho_c, A, psi):
    """
    A à employer dans la forme normalisée (rho_i,c = 0).

    En posant y = (rho_i - rho_i_c)/(rho_i_v - rho_i_c) dans l'équation native
        d(rho_i)/dt = -A · rho_i_v · [(rho_i - rho_i_c)/rho_i_v]^psi · exp(...)
    il vient
        dy/dt = -A · [(rho_i_v - rho_i_c)/rho_i_v]^(psi-1) · y^psi · exp(...)
    Le crochet vaut 1 si rho_i_c = 0 — sinon il faut le porter dans A.
    """
    return A * ((rho_v - rho_c) / rho_v) ** (psi - 1.0)


# ---------------------------------------------------------------------------
# Émission
# ---------------------------------------------------------------------------

def show_native(mat):
    print("=" * 79)
    print("  JEU 0 — FORME NATIVE  (rho_i,v ; rho_i,c)  — telle que publiée")
    print("=" * 79)
    print("""
      d(rho_i)/dt = - A_i · rho_i,v · [(rho_i - rho_i,c)/rho_i,v]^psi_i
                            · exp(-E_i/(R·T))       pour T > T_reac, 0 sinon
""")
    print("  (a) base RÉSINE — kg/m³ de matrice, ce que donne le classeur :\n")
    print(f"      {'phase':>6} {'rho_i,v':>9} {'rho_i,c':>9} {'A [1/s]':>12} "
          f"{'E/R [K]':>13} {'psi':>5} {'T_reac [K]':>11}")
    for nm, rv, rc, A, er, ps, tr in GOLDSTEIN:
        print(f"      {nm:>6} {rv:9.1f} {rc:9.1f} {A:12.4e} {er:13.6f} "
              f"{ps:5.1f} {tr:11.4f}")
    print(f"      {'C':>6} {RHO_FIBER:9.1f} {RHO_FIBER:9.1f} {0.0:12.4e} "
          f"{'—':>13} {'—':>5} {'—':>11}   fibres, inertes")
    print(f"\n      somme rho_i,v = {sum(g[1] for g in GOLDSTEIN):.0f} kg/m³ "
          f"= densité intrinsèque de la résine")

    e = mat['eps_m']
    print(f"\n  (b) base COMPOSITE — multiplier par eps_m = {e:.6g} :\n")
    print(f"      {'phase':>6} {'rho_i,v':>9} {'rho_i,c':>9} {'A [1/s]':>12} "
          f"{'E/R [K]':>13} {'psi':>5} {'T_reac [K]':>11}")
    for nm, rv, rc, A, er, ps, tr in GOLDSTEIN:
        print(f"      {nm:>6} {e*rv:9.1f} {e*rc:9.1f} {A:12.4e} {er:13.6f} "
              f"{ps:5.1f} {tr:11.4f}")
    print(f"      {'C':>6} {mat['eps_f']*RHO_FIBER:9.1f} {mat['eps_f']*RHO_FIBER:9.1f} "
          f"{0.0:12.4e} {'—':>13} {'—':>5} {'—':>11}   fibres")
    tv = sum(e*g[1] for g in GOLDSTEIN) + mat['eps_f']*RHO_FIBER
    tc = sum(e*g[2] for g in GOLDSTEIN) + mat['eps_f']*RHO_FIBER
    print(f"\n      somme vierge {tv:.1f} (attendu {mat['rho_v']:.0f})   "
          f"somme char {tc:.1f} (attendu {mat['rho_c']:.0f})")
    print("\n      Les A sont INCHANGÉS : l'équation native est invariante par")
    print("      mise à l'échelle commune de (rho_i, rho_i,v, rho_i,c).")


def show_correction():
    print("\n" + "=" * 79)
    print("  LA CORRECTION QU'EXIGENT LES JEUX A ET B")
    print("=" * 79)
    print("""
  Les formes « en f » et « en Delta_rho » posent implicitement rho_i,c = 0.
  Or la phase B du TACOT a rho_c = 600. Le facteur à reporter dans A vaut
  [(rho_i,v - rho_i,c)/rho_i,v]^(psi-1) :
""")
    print(f"      {'phase':>6} {'(rho_v-rho_c)/rho_v':>21} {'facteur':>10} "
          f"{'A publié':>12} {'A à employer':>14}")
    out = []
    for nm, rv, rc, A, er, ps, tr in GOLDSTEIN:
        fac = ((rv - rc) / rv) ** (ps - 1.0)
        Ac = corrected_A(rv, rc, A, ps)
        out.append(Ac)
        print(f"      {nm:>6} {(rv-rc)/rv:21.6f} {fac:10.6f} {A:12.4e} "
              f"{Ac:14.4e}")
    print("\n      La phase A est intacte (son rho_c était déjà nul) ; seule la")
    print("      phase B change, d'un facteur 9. Oublier la correction rend la")
    print("      pyrolyse 9 fois trop rapide sur cette réaction — voir la")
    print("      vérification plus bas.")
    return out


def show_deck_a(mat, Ac):
    print("\n" + "=" * 79)
    print("  JEU A — ARRHENIUS + F_i   (fractions, base composite)")
    print("=" * 79)
    print("""
      dx_i/dt = - A_i · x_i^psi_i · exp(-E_i/(R·T))     x_i(0) = 1
      rho(T)  = rho_vierge · [ 1 - somme_i F_i · (1 - x_i) ]
""")
    print(f"      rho_vierge      {mat['rho_v']:9.2f} kg/m³")
    print(f"      rho_char        {mat['rho_c']:9.2f} kg/m³")
    print(f"      porosité vierge {1-mat['eps_f']-mat['eps_m']:9.4f}")
    print()
    print(f"      {'i':>2} {'F_i [-]':>12} {'A_i [1/s]':>12} {'E/R [K]':>13} "
          f"{'psi':>5} {'T_reac [K]':>11}")
    Fs = []
    for (nm, rv, rc, A, er, ps, tr), a in zip(GOLDSTEIN, Ac):
        F = mat['eps_m'] * (rv - rc) / mat['rho_v']
        Fs.append(F)
        print(f"      {nm:>2} {F:12.9f} {a:12.4e} {er:13.6f} {ps:5.1f} "
              f"{tr:11.4f}")
    print(f"      {'':>2} {sum(Fs):12.9f}   <- somme")
    print(f"      contrôle : 1 - {sum(Fs):.9f} = {1-sum(Fs):.9f}   "
          f"rho_c/rho_v = {mat['rho_c']/mat['rho_v']:.9f}")
    return Fs


def show_deck_b(mat, Ac):
    print("\n" + "=" * 79)
    print("  JEU B — ARRHENIUS + Delta_rho_i   (kg/m³ de composite)")
    print("=" * 79)
    print("""
      d(rho_i)/dt = - A_i · Drho_i · (rho_i/Drho_i)^psi_i · exp(-E_i/(R·T))
      rho_i(0) = Drho_i        rho(T) = rho_inerte + somme_i rho_i
""")
    print(f"      rho_vierge         {mat['rho_v']:9.2f} kg/m³")
    print(f"      rho_inerte (=char) {mat['rho_c']:9.2f} kg/m³"
          f"   dont fibres {mat['eps_f']*RHO_FIBER:.0f} et char de résine "
          f"{mat['rho_c'] - mat['eps_f']*RHO_FIBER:.0f}")
    print()
    print(f"      {'i':>2} {'Drho_i [kg/m3]':>15} {'A_i [1/s]':>12} "
          f"{'E/R [K]':>13} {'psi':>5} {'T_reac [K]':>11}")
    Ds = []
    for (nm, rv, rc, A, er, ps, tr), a in zip(GOLDSTEIN, Ac):
        D = mat['eps_m'] * (rv - rc)
        Ds.append(D)
        print(f"      {nm:>2} {D:15.6f} {a:12.4e} {er:13.6f} {ps:5.1f} "
              f"{tr:11.4f}")
    print(f"      {'':>2} {sum(Ds):15.6f}   <- somme = gaz total libéré")
    print(f"      contrôle : {mat['rho_c']:.2f} + {sum(Ds):.2f} = {mat['rho_c']+sum(Ds):.2f}"
          f"   (rho_vierge = {mat['rho_v']:.2f})")
    print("\n      Si le code réclame un couple (rho_i,v ; rho_i,c) : soit vous")
    print("      employez le JEU 0 tel quel avec les A publiés, soit vous posez")
    print("      rho_i,v = Drho_i et rho_i,c = 0 AVEC les A corrigés ci-dessus.")
    print("      Panacher les deux est l'erreur classique.")
    return Ds


# ---------------------------------------------------------------------------
# Vérification
# ---------------------------------------------------------------------------

def _grid(beta, t0, t1, n):
    b = beta / 60.0
    dT = (t1 - t0) / n
    return dT, dT / b


def integrate_native(mat, beta, t0=300.0, t1=1400.0, n=200000):
    dT, dt = _grid(beta, t0, t1, n)
    r = [mat['eps_m'] * g[1] for g in GOLDSTEIN]
    inert = mat['eps_f'] * RHO_FIBER
    Ts, rho = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        rho.append(inert + sum(r))
        for j, (nm, rv, rc, A, er, ps, tr) in enumerate(GOLDSTEIN):
            rvj, rcj = mat['eps_m'] * rv, mat['eps_m'] * rc
            if T > tr and r[j] > rcj:
                r[j] = max(rcj, r[j] - A * rvj * ((r[j] - rcj) / rvj) ** ps
                           * math.exp(-er / T) * dt)
    return Ts, rho


def integrate_a(mat, beta, Fs, Ac, t0=300.0, t1=1400.0, n=200000):
    dT, dt = _grid(beta, t0, t1, n)
    x = [1.0] * len(GOLDSTEIN)
    Ts, rho = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        rho.append(mat['rho_v'] * (1.0 - sum(F * (1.0 - xi) for F, xi in zip(Fs, x))))
        for j, (nm, rv, rc, A, er, ps, tr) in enumerate(GOLDSTEIN):
            if T > tr and x[j] > 0.0:
                x[j] = max(0.0, x[j] - Ac[j] * x[j] ** ps
                           * math.exp(-er / T) * dt)
    return Ts, rho


def integrate_b(mat, beta, Ds, Ac, t0=300.0, t1=1400.0, n=200000):
    dT, dt = _grid(beta, t0, t1, n)
    r = list(Ds)
    inert = mat['rho_v'] - sum(Ds)
    Ts, rho = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        rho.append(inert + sum(r))
        for j, (nm, rv, rc, A, er, ps, tr) in enumerate(GOLDSTEIN):
            if T > tr and r[j] > 0.0:
                r[j] = max(0.0, r[j] - Ac[j] * Ds[j] * (r[j] / Ds[j]) ** ps
                           * math.exp(-er / T) * dt)
    return Ts, rho


def verify(mat, beta, Fs, Ds, Ac):
    print("\n" + "=" * 79)
    print(f"  VÉRIFICATION — les trois jeux décrivent-ils le même matériau ?")
    print("=" * 79)
    print(f"  Intégrations indépendantes, ATG à {beta:.0f} K/min.\n")
    Tn, Rn = integrate_native(mat, beta)
    _, Ra = integrate_a(mat, beta, Fs, Ac)
    _, Rb = integrate_b(mat, beta, Ds, Ac)
    A_nc = [g[3] for g in GOLDSTEIN]                     # A NON corrigés
    _, Rx = integrate_b(mat, beta, Ds, A_nc)

    print(f"  {'T [°C]':>7} | {'JEU 0 natif':>12} {'JEU A':>11} {'JEU B':>11}"
          f" | {'sans correction':>16}")
    for Tc in (300, 400, 500, 550, 600, 700, 900):
        i = min(range(len(Tn)), key=lambda k: abs(Tn[k] - (Tc + 273.15)))
        print(f"  {Tc:>7} | {Rn[i]:12.4f} {Ra[i]:11.4f} {Rb[i]:11.4f}"
              f" | {Rx[i]:16.4f}")
    da = max(abs(a - b) for a, b in zip(Rn, Ra))
    db = max(abs(a - b) for a, b in zip(Rn, Rb))
    dx = max(abs(a - b) for a, b in zip(Rn, Rx))
    print(f"\n  écart max au JEU 0 :  JEU A {da:.2e}   JEU B {db:.2e} kg/m³")
    ok = da < 1e-8 and db < 1e-8
    print(f"  -> {'LES TROIS JEUX SONT ÉQUIVALENTS' if ok else 'ÉCHEC'}")
    print(f"\n  Et si l'on oublie la correction sur A_B : écart {dx:.1f} kg/m³,")
    print(f"  soit {100*dx/(mat['rho_v']-mat['rho_c']):.0f} % du gaz total libéré. L'erreur est")
    print("  du premier ordre — ce n'est pas un raffinement.")
    return ok


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("material", nargs="?", default="TACOT",
                    choices=sorted(MATERIALS), help="matériau à traiter")
    ap.add_argument("--beta", type=float, default=20.0,
                    help="vitesse de chauffe de l'ATG simulée [K/min]")
    args = ap.parse_args()
    mat = microstructure(args.material)

    print()
    print(f"  MISES EN DONNÉES DE PYROLYSE — {mat['name']}")
    print(f"  {mat['note']}")
    print("  Cinétique de Goldstein via TACOT_3.0.xls, 'Pyrolysis model'!A11:G20")
    if mat["name"] != "TACOT":
        print("  (mêmes fibres et même résine que le TACOT : seul eps_m change)")
    print()
    show_native(mat)
    Ac = show_correction()
    Fs = show_deck_a(mat, Ac)
    Ds = show_deck_b(mat, Ac)
    ok = verify(mat, args.beta, Fs, Ds, Ac)

    print("\n" + "=" * 79)
    print("  RAPPEL — ce qui accompagne le bloc cinétique")
    print("=" * 79)
    print("      gaz de pyrolyse   " + ", ".join(f"{e}:{v}" for e, v in PYRO_GAS.items())
          + "   (fractions molaires)")
    print("      char              C:1.0")
    print("      bord de couche l. N:0.79, O:0.21  (air)")
    print(f"      rho_v / rho_c     {mat['rho_v']:.2f} / {mat['rho_c']:.2f} kg/m³")
    print(f"      eps_f / eps_m     {mat['eps_f']:.4f} / {mat['eps_m']:.4f}"
          f"   porosité {1-mat['eps_f']-mat['eps_m']:.4f}")
    k = (mat['rho_v'] - mat['rho_c']) / mat['rho_c']
    print(f"      couplage k        {k:.4f}   = B'g/B'c en ablation stationnaire")
    print("\n      Le classeur d'origine ne porte AUCUNE enthalpie de pyrolyse ;")
    print("      les -4e6 J/kg qu'on voit parfois viennent de la base VKI.")
    print("=" * 79)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
