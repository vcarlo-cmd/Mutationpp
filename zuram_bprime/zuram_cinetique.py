#!/usr/bin/env python
"""
Paramètres d'Arrhenius du ZURAM : peut-on les transposer aux variantes XX/YY ?

Réponse courte : A, E/R et m ne bougent PAS, et f non plus. Ce sont des
propriétés de la RÉSINE, or la résine est la même. Seul change le RACCORD au
composite — la quantité de matière qui se décompose, pas sa manière de le faire.

Ce script établit trois choses :

  1. LA FORME NORMALISÉE. Le classeur donne f, log10A, E/R, m mais AUCUN
     rho_i_v. On montre qu'une seule lecture est auto-cohérente : chaque
     réaction consomme intégralement sa pseudo-phase (rho_i_c = 0), le char
     étant porté par une phase inerte séparée. Sous cette lecture le système
     est invariant d'échelle et (f, A, E/R, m) suffit.

  2. LA TRANSPOSITION. Tables prêtes à recopier, dans les deux conventions
     usuelles : base RÉSINE (inchangée) et base COMPOSITE (ré-échelonnée).

  3. LA VÉRIFICATION. Simulation d'ATG à vitesse de chauffe constante pour
     chaque variante : les courbes normalisées par la masse de résine se
     superposent exactement, les pics de DTG sont aux mêmes températures, et
     la perte du composite vaut bien w_résine x 0.380234.

  4. UNE CONTRE-ÉPREUVE. Le classeur réécrit aussi la cinétique du TACOT dans
     cette forme en f, et là on dispose de l'original : la transcription
     confirme que les f sont bien des grandeurs de RÉSINE, mais révèle que la
     lecture naïve (rho_i,c = 0) ne reproduit pas la cinétique de départ.

Source : 'ZURAM_official'!A73:F83 du classeur
Tacot_Zuram_Calcarb_database_v4.3.1.ods. Voir `variantes_zuram.md`.

Usage :
    python zuram_cinetique.py
    python zuram_cinetique.py 14/40 18/80 --beta 10
"""

import argparse
import math

from zuram_variantes import KINETICS, RHO_RESIN, CHAR_YIELD, parse, variant

SUM_F = sum(f for f, _, _, _ in KINETICS)


# ---------------------------------------------------------------------------
# 1. Forme normalisée
# ---------------------------------------------------------------------------

def show_form():
    print("=" * 79)
    print("  1. QUELLE ÉQUATION, AU JUSTE ?")
    print("=" * 79)
    print("""
  Le classeur écrit ses paramètres dans la forme dite « TACOT-like »
  (note de version 1.4.0), celle de la feuille 'Pyrolysis model' du TACOT :

      d(rho_i)/dt = - A_i · rho_i_v · [(rho_i - rho_i_c)/rho_i_v]^m
                             · exp(-E_i/(R·T))

  Mais il ne fournit QUE f_i, A_i, E_i/R et m_i — jamais rho_i_v ni rho_i_c
  séparément. Or l'équation ci-dessus n'est pas invariante d'échelle : en
  posant y = (rho_i - rho_i_c)/(rho_i_v - rho_i_c), on obtient

      dy/dt = - A_i · [(rho_i_v - rho_i_c)/rho_i_v]^(m-1) · y^m · exp(-E/RT)

  Le facteur entre crochets ne disparaît QUE si rho_i_c = 0.

  -> Une seule lecture rend la donnée du classeur complète : chaque réaction
     consomme INTÉGRALEMENT sa pseudo-phase (rho_i_c = 0), le char résiduel
     étant porté par une phase inerte distincte. On a alors, avec
     x_i = rho_i/rho_i_v dans [0, 1] :

         dx_i/dt = - A_i · x_i^m · exp(-E_i/(R·T))          (1)

     rho_résine(t) = rho_résine,v · [ (1 - somme f_i) + somme f_i · x_i ]

  C'est la forme employée ci-dessous. Si votre code attend un rho_i_c non nul,
  il faut corriger A par le facteur [(rho_i_v - rho_i_c)/rho_i_v]^(m-1), sans
  quoi les vitesses sont fausses — avec m entre 2.6 et 4.6, l'erreur est
  d'ordre un.
""")
    print(f"  Contrôle : somme des f = {SUM_F:.6f}, "
          f"donc rendement en char {1-SUM_F:.6f}")
    print(f"  'ZURAM_official' donne {CHAR_YIELD:.6f} — écart "
          f"{abs(1-SUM_F-CHAR_YIELD):.1e}")


# ---------------------------------------------------------------------------
# 2. Transposition
# ---------------------------------------------------------------------------

def show_tables(variants, labels):
    print("\n" + "=" * 79)
    print("  2. CE QU'IL FAUT ÉCRIRE POUR CHAQUE VARIANTE")
    print("=" * 79)

    print("\n  (a) BASE RÉSINE — la convention du classeur. RIEN NE CHANGE.")
    print(f"      Identique pour 14/40, 18/50, 18/80 et toute autre variante :\n")
    print(f"      {'réaction':>9} {'f [-]':>12} {'log10 A':>9} {'A [1/s]':>12}"
          f" {'E/R [K]':>11} {'m':>6}")
    for i, (f, la, er, m) in enumerate(KINETICS, 1):
        print(f"      {i:>9} {f:12.6f} {la:9.2f} {10**la:12.3e} {er:11.1f} {m:6.2f}")
    print("\n      A, E/R et m décrivent la CHIMIE du novolac ; f décrit la")
    print("      part de résine que chaque réaction consomme. Ni l'une ni")
    print("      l'autre ne connaît la préforme.")

    print("\n  (b) BASE COMPOSITE — si le code raisonne sur la densité totale.")
    print("      Seules les fractions de perte se ré-échelonnent : F_i = f_i · w.\n")
    print(f"      {'variante':>10} {'w_résine':>10}"
          + "".join(f"{'F' + str(i):>12}" for i in range(1, 5)) + f"{'somme':>11}")
    for lab, v in zip(labels, variants):
        Fs = [f * v["w_resin"] for f, _, _, _ in KINETICS]
        print(f"      {lab:>10} {v['w_resin']:10.3f}"
              + "".join(f"{F:12.6f}" for F in Fs) + f"{sum(Fs):11.6f}")

    print("\n  (c) BASE VOLUMIQUE — si le code veut des rho_i_v en kg/m³ de")
    print("      COMPOSITE (style TACOT 'Pyrolysis model'!B18:C20).")
    print("      rho_i_v = f_i · eps_m · rho_résine,intr,  rho_i_c = 0\n")
    print(f"      {'variante':>10} {'eps_m':>8}"
          + "".join(f"{'rho' + str(i) + '_v':>11}" for i in range(1, 5))
          + f"{'char':>11}{'fibres':>11}")
    for lab, v in zip(labels, variants):
        rv = [f * v["eps_r"] * RHO_RESIN for f, _, _, _ in KINETICS]
        char = (1.0 - SUM_F) * v["eps_r"] * RHO_RESIN
        print(f"      {lab:>10} {v['eps_r']:8.4f}"
              + "".join(f"{r:11.2f}" for r in rv)
              + f"{char:11.2f}{v['m_f']:11.2f}")
    print("\n      Contrôle : somme(rho_i_v) + char + fibres = rho vierge")
    for lab, v in zip(labels, variants):
        tot = SUM_F * v["eps_r"] * RHO_RESIN + (1 - SUM_F) * v["eps_r"] * RHO_RESIN \
              + v["m_f"]
        print(f"      {lab:>10} -> {tot:8.2f}   (attendu {v['rho_v']:8.2f})")

    print("\n      ATTENTION à la convention de pyrolyse. Le classeur suppose une")
    print("      densité intrinsèque de matrice INCHANGÉE avec perte de VOLUME")
    print("      (anomalie ouverte n° 3), alors que la cinétique décrit une perte")
    print("      de DENSITÉ à volume constant. Les deux donnent le même rho_char")
    print("      composite — mais ne pas mélanger les deux dans un même calcul.")


# ---------------------------------------------------------------------------
# 3. Vérification par simulation d'ATG
# ---------------------------------------------------------------------------

def tga_resin(beta, t_start=300.0, t_end=1400.0, n=200000):
    """
    Base RÉSINE : intègre dx_i/dt = -A_i x_i^m exp(-E_i/RT) à vitesse de
    chauffe constante beta [K/min]. Retourne (T, perte, vitesse de perte).
    """
    b = beta / 60.0
    dT = (t_end - t_start) / n
    dt = dT / b
    x = [1.0] * len(KINETICS)
    Ts, loss, rate = [], [], []
    prev = 0.0
    for k in range(n + 1):
        T = t_start + k * dT
        L = sum(f * (1.0 - xi) for (f, _, _, _), xi in zip(KINETICS, x))
        Ts.append(T)
        loss.append(L)
        rate.append((L - prev) / dT if k else 0.0)
        prev = L
        for j, (f, la, er, m) in enumerate(KINETICS):
            if x[j] > 0.0:
                x[j] = max(0.0, x[j] - 10.0 ** la * x[j] ** m
                           * math.exp(-er / T) * dt)
    return Ts, loss, rate


def tga_composite(v, beta, t_start=300.0, t_end=1400.0, n=200000):
    """
    Base COMPOSITE : intègre la MÊME équation mais sur les rho_i exprimés en
    kg/m³ de composite, tels que la § 2(c) les prescrit —

        rho_i_v = f_i · eps_m · rho_résine,intr ,   rho_i_c = 0
        d(rho_i)/dt = -A_i · rho_i_v · (rho_i/rho_i_v)^m · exp(-E_i/RT)

    Retourne (T, rho(T)). Aucune référence à la courbe résine : c'est ce qui
    rend la comparaison non circulaire.
    """
    b = beta / 60.0
    dT = (t_end - t_start) / n
    dt = dT / b
    rho_iv = [f * v["eps_r"] * RHO_RESIN for f, _, _, _ in KINETICS]
    rho = list(rho_iv)
    inert = v["m_f"] + (1.0 - SUM_F) * v["eps_r"] * RHO_RESIN
    Ts, rhos = [], []
    for k in range(n + 1):
        T = t_start + k * dT
        Ts.append(T)
        rhos.append(inert + sum(rho))
        for j, (f, la, er, m) in enumerate(KINETICS):
            if rho[j] > 0.0:
                d = 10.0 ** la * rho_iv[j] * (rho[j] / rho_iv[j]) ** m \
                    * math.exp(-er / T) * dt
                rho[j] = max(0.0, rho[j] - d)
    return Ts, rhos


def show_verification(variants, labels, beta):
    print("\n" + "=" * 79)
    print(f"  3. VÉRIFICATION — ATG simulée à {beta:.0f} K/min")
    print("=" * 79)
    Ts, loss, rate = tga_resin(beta)

    print("\n  On intègre SÉPARÉMENT la base résine et, pour chaque variante,")
    print("  la base composite avec ses propres rho_i_v (§ 2c). Les deux")
    print("  doivent coïncider une fois normalisées — sinon la transposition")
    print("  est fausse.\n")

    comp = [tga_composite(v, beta) for v in variants]

    print(f"  {'T [°C]':>7} | {'perte résine':>13} |"
          + "".join(f"{'  ' + lab:>16}" for lab in labels))
    print(f"  {'':7} | {'':13} |"
          + "".join(f"{'composite':>10}{'attendu':>6}" for _ in labels))
    for Tc in (400, 500, 600, 650, 700, 800, 1000):
        i = min(range(len(Ts)), key=lambda k: abs(Ts[k] - (Tc + 273.15)))
        cells = ""
        for v, (_, rhos) in zip(variants, comp):
            got = (v["rho_v"] - rhos[i]) / v["rho_v"]
            cells += f"{100*got:9.3f}%{100*loss[i]*v['w_resin']:6.2f}"
        print(f"  {Tc:>7} | {100*loss[i]:12.3f} % |" + cells)

    worst = 0.0
    for v, (_, rhos) in zip(variants, comp):
        for i in range(0, len(Ts), 97):
            got = (v["rho_v"] - rhos[i]) / v["rho_v"]
            worst = max(worst, abs(got - loss[i] * v["w_resin"]))
    ok = worst < 1e-9
    print(f"\n  écart maximal entre les deux intégrations : {worst:.2e}")
    print(f"  -> {'TRANSPOSITION EXACTE' if ok else 'ÉCHEC'} : "
          f"(rho_v - rho)/rho_v = w_résine · perte_résine(T)")

    imax = max(range(1, len(rate)), key=lambda k: rate[k])
    print(f"\n  Pic de DTG à T = {Ts[imax]:.1f} K ({Ts[imax]-273.15:.1f} °C),")
    print("  identique pour toutes les variantes : c'est la même résine.")
    print("  Sur le composite, ce pic est simplement moins marqué — les fibres")
    print("  Calcarb ne pyrolysent pas et diluent le signal d'un facteur w.")

    print(f"\n  QUEUE ALGÉBRIQUE. À 1400 K la perte de la résine atteint")
    print(f"  {100*loss[-1]:.3f} % et non les {100*SUM_F:.3f} % de la somme des f.")
    print("  Ce n'est pas une erreur d'intégration (le résultat est convergé en")
    print("  pas) mais une propriété du modèle : avec m > 1, dx/dt = -A·x^m")
    print("  donne x ~ t^(-1/(m-1)), qui ne s'annule jamais. Le rendement en")
    print("  char de 0.6198 est une valeur ASYMPTOTIQUE ; à 1400 K le modèle")
    print("  en rend 0.6244. À prendre en compte si l'on compare à une ATG")
    print("  qui, elle, s'arrête à une température finie.")
    return ok


# ---------------------------------------------------------------------------
# 4. Contre-épreuve sur le TACOT : la forme en f est-elle fidèle ?
# ---------------------------------------------------------------------------

# TACOT_3.0.xls, 'Pyrolysis model'!B18:G19 — cinétique de Goldstein d'origine
TACOT_GOLDSTEIN = [   # rho_i_v, rho_i_c, A, E/R, T_reac
    (300.0,   0.0, 1.2e4,  8555.555555555555, 333.3333333333333),
    (900.0, 600.0, 4.48e9, 20444.444444444445, 555.5555555555555),
]
TACOT_RHO_RESIN = 1200.0
# Tacot_Zuram_Calcarb_database, onglet 'TACOT_3_0'!A110:F111 — même cinétique
# réécrite dans la forme en f, avec les A repris VERBATIM
TACOT_VKI_F = [  # f, A, E/R, T_reac
    (0.25, 1.2e4,  8555.555555555555, 333.3333333333333),
    (0.25, 4.48e9, 20444.444444444445, 555.5555555555555),
]


def _tga_goldstein(beta, t0=300.0, t1=1400.0, n=200000):
    b, dT = beta / 60.0, (t1 - t0) / n
    dt = dT / b
    r = [p[0] for p in TACOT_GOLDSTEIN]
    Ts, loss = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        loss.append((TACOT_RHO_RESIN - sum(r)) / TACOT_RHO_RESIN)
        for j, (rv, rc, A, er, tr) in enumerate(TACOT_GOLDSTEIN):
            if T > tr and r[j] > rc:
                r[j] = max(rc, r[j] - A * rv * ((r[j] - rc) / rv) ** 3
                           * math.exp(-er / T) * dt)
    return Ts, loss


def _tga_fform(beta, scales=None, t0=300.0, t1=1400.0, n=200000):
    """scales : facteur correctif par réaction (1.0 = A du classeur tel quel)."""
    if scales is None:
        scales = [1.0] * len(TACOT_VKI_F)
    b, dT = beta / 60.0, (t1 - t0) / n
    dt = dT / b
    x = [1.0] * len(TACOT_VKI_F)
    Ts, loss = [], []
    for k in range(n + 1):
        T = t0 + k * dT
        Ts.append(T)
        loss.append(sum(f * (1.0 - xi) for (f, _, _, _), xi in zip(TACOT_VKI_F, x)))
        for j, (f, A, er, tr) in enumerate(TACOT_VKI_F):
            if T > tr and x[j] > 0.0:
                x[j] = max(0.0, x[j] - A * scales[j] * x[j] ** 3
                           * math.exp(-er / T) * dt)
    return Ts, loss


def _peak(Ts, loss):
    d = [loss[k + 1] - loss[k] for k in range(len(loss) - 1)]
    return Ts[max(range(len(d)), key=lambda k: d[k])] - 273.15


def check_tacot_transcription(beta):
    print("\n" + "=" * 79)
    print("  4. CONTRE-ÉPREUVE — la forme en f est-elle fidèle ? (cas du TACOT)")
    print("=" * 79)
    print("""
  Le classeur VKI contient aussi un onglet 'TACOT_3_0' où la cinétique de
  Goldstein a été réécrite dans la MÊME forme en f que le ZURAM. Comme on
  dispose ici de l'original, on peut vérifier la transcription.

  Original ('Pyrolysis model'!B18:C19) :
      phase A : rho_v = 300, rho_c =   0,  A = 1.2e4
      phase B : rho_v = 900, rho_c = 600,  A = 4.48e9
  VKI ('TACOT_3_0'!A110:C111) :
      f1 = 0.25, A = 1.2e4      f2 = 0.25, A = 4.48e9      (A repris verbatim)
""")
    print(f"  f1 = (300-0)/1200   = {(300-0)/1200}")
    print(f"  f2 = (900-600)/1200 = {(900-600)/1200}")
    print("  -> les f sont bien normalisés par la densité INTRINSÈQUE DE LA")
    print("     RÉSINE (1200 kg/m3), et leur somme 0.50 est son rendement en char.")
    print("     Sur ce point, la transcription est juste : ce sont des grandeurs")
    print("     de RÉSINE, jamais de composite.\n")

    print("  Mais la phase B avait un rho_c NON NUL, que la forme en f ne sait")
    print("  pas exprimer. Vitesse initiale, hors exponentielle :\n")
    for lab, (rv, rc, A, _, _) in zip(("phase A", "phase B"), TACOT_GOLDSTEIN):
        orig, ff = A * rv * ((rv - rc) / rv) ** 3, A * (rv - rc)
        print(f"      {lab} : original {orig:11.4e}   forme f {ff:11.4e}"
              f"   rapport {ff/orig:5.2f}")
    print(f"\n  facteur attendu [(rho_v-rho_c)/rho_v]^(1-m) = (300/900)^-2 = 9")
    print(f"  A2 qu'il aurait fallu écrire : 4.48e9/9 = {4.48e9/9:.4e}\n")

    Tg, Lg = _tga_goldstein(beta)
    Tv, Lv = _tga_fform(beta)                        # A du classeur, verbatim
    Tc, Lc = _tga_fform(beta, [1.0, 1.0 / 9.0])      # seule la phase B corrigée
    print(f"  ATG de la résine du TACOT à {beta:.0f} K/min :\n")
    print(f"  {'T [°C]':>7} | {'Goldstein':>11} {'VKI verbatim':>14} {'VKI A2/9':>11}")
    for T in (300, 400, 500, 550, 600, 700, 900):
        i = min(range(len(Tg)), key=lambda k: abs(Tg[k] - (T + 273.15)))
        print(f"  {T:>7} | {100*Lg[i]:10.2f}% {100*Lv[i]:13.2f}% {100*Lc[i]:10.2f}%")
    dv = max(abs(a - b) for a, b in zip(Lg, Lv))
    dc = max(abs(a - b) for a, b in zip(Lg, Lc))
    print(f"\n  écart max / Goldstein : verbatim {100*dv:5.2f} points, "
          f"avec A2/9 {100*dc:5.2f} points")
    print(f"  pic de DTG : Goldstein {_peak(Tg, Lg):.1f} °C | "
          f"VKI verbatim {_peak(Tv, Lv):.1f} °C | VKI A2/9 {_peak(Tc, Lc):.1f} °C")
    print("""
  CONCLUSION. Lue naïvement (rho_i,c = 0), la forme en f du classeur ne
  reproduit PAS la cinétique du TACOT : 9.9 points d'écart et un pic de DTG
  décalé de 63 °C. Deux lectures possibles, et le classeur ne tranche pas :
     (a) il manque un rho_i,c, que f seul ne peut pas porter ;
     (b) ou A2 aurait dû être divisé par 9 à la transcription.

  POUR LE ZURAM, on ne peut pas trancher de la même façon : ses f viennent
  d'un ajustement direct sur ATG (source [5] du classeur), pas de la
  conversion d'un jeu (rho_v, rho_c) préexistant. La lecture normalisée est
  la seule que la donnée fournie rende auto-suffisante, et elle place le pic
  de DTG à 557 °C — cohérent avec le maximum de vitesse « ~600 °C » rapporté
  par Torres-Herrador. Mais c'est un faisceau d'indices, pas une preuve :
  à confirmer auprès du VKI avant tout calcul de production.
""")
    return dc < 1e-6


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*", default=["14/40", "18/50", "18/80"])
    ap.add_argument("--beta", type=float, default=20.0,
                    help="vitesse de chauffe de l'ATG simulée [K/min]")
    args = ap.parse_args()

    labels = args.variants
    variants = [variant(*parse(a)) for a in labels]

    print()
    print("  PARAMÈTRES D'ARRHENIUS DU ZURAM — TRANSPOSITION AUX VARIANTES XX/YY")
    print()
    show_form()
    show_tables(variants, labels)
    ok = show_verification(variants, labels, args.beta)
    ok &= check_tacot_transcription(args.beta)

    print("\n" + "=" * 79)
    print("  CONCLUSION")
    print("=" * 79)
    print("  Il n'y a RIEN à extrapoler sur A, E/R et m : ce sont des grandeurs")
    print("  de la chimie du novolac, pas du composite. f non plus, dès lors")
    print("  qu'on l'exprime sur la résine — ce que le classeur fait déjà")
    print("  (note B77). La seule opération est le raccord au composite :")
    print("      F_i = f_i · YY/100          (base massique)")
    print("      rho_i_v = f_i · eps_m · rho_résine,intr   (base volumique)")
    print()
    print("  Cette invariance suppose une résine identiquement RÉTICULÉE. C'est")
    print("  l'hypothèse fragile pour 18/80 : voir `variantes_zuram.md` § 4.")
    print("=" * 79)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
