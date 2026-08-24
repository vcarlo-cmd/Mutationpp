#!/usr/bin/env python
"""
Variantes ZURAM XX/YY : ce que le nom impose, et ce qu'il ne touche pas.

En lisant le nom « ZURAM 18/50 » comme (cf. `resine_zuram.md` § 1) :

    XX = masse volumique apparente de la préforme, en 0.01 g/cm³
    YY = teneur en résine du composite vierge, en % MASSE

tout le bilan de phase solide se déduit de deux nombres, à fibres et résine
inchangées. Ce script explicite ce qui bouge et ce qui ne bouge pas.

Le résultat central est une SÉPARATION DES RÔLES :

    XX  ne fixe QUE l'échelle (densités, masse) — un simple facteur
    YY  fixe TOUTE la thermochimie de la réponse matériau (couplage k,
        perte de masse en ATG, porosité)

et, dans les deux cas, la TABLE B' EST INCHANGÉE : elle ne consomme que des
compositions élémentaires normalisées, or ni la résine ni les fibres ne
changent. C'est le même argument que `../tacot_bprime/cph70_vs_tacot.md`,
vérifié là-bas numériquement à 0.000e+00 près.

Usage :
    python zuram_variantes.py
    python zuram_variantes.py 14/40 18/80 20/65
"""

import sys

# --- constantes matériau, de 'ZURAM_official' (cf. resine_zuram.md) ---------
RHO_FIBER = 1577.2413793      # kg/m³, densité intrinsèque des fibres Calcarb
RHO_RESIN = 1315.0719638      # kg/m³, densité intrinsèque de la résine vierge
# Rendement en char, en PLEINE PRÉCISION depuis les fractions volumiques du
# classeur ('ZURAM_official'!F12:F13). La seconde voie du classeur — 1 moins la
# somme des f de la cinétique — donne la même valeur à 3.3e-16 près. Un arrondi
# à 6 chiffres suffit à faire diverger rho_char et rho_vierge·(1-somme F) au
# 7e chiffre : ne pas tronquer.
CHAR_YIELD = 0.100853796331154 / 0.1627287372   # = 0.619766355141076
# hypothèse du classeur : densité intrinsèque de matrice INCHANGÉE à la
# pyrolyse (pure perte de volume). C'est l'anomalie ouverte n° 3.
RHO_RESIN_CHAR = RHO_RESIN

# gaz de pyrolyse et char : INDÉPENDANTS de XX et YY
PYRO_GAS = {"C": 0.171, "H": 0.722, "O": 0.107}
CHAR = {"C": 1.0}

# Cinétique de pyrolyse, 'ZURAM_official'!A79:F83.
# Les f sont DÉJÀ ramenées à la résine seule : la note B77 précise que les F
# de la source, exprimées sur le composite, ont été divisées par la fraction
# massique de résine. Elles sont donc transposables telles quelles à toute
# variante XX/YY — c'est le même polymère.
KINETICS = [  # (f, log10 A, E/R [K], m)
    (0.035070093457862, 5.33, 8178.52014290537, 4.30),
    (0.0276869158877858, 8.69, 16068.3866337082, 3.70),
    (0.0959813084109907, 10.60, 21612.9422011779, 2.57),
    (0.221495327102286, 11.67, 26423.8364028869, 4.63),
]

# référence mesurée (et non nominale) du ZURAM_official
REF = dict(rho_v=418.998084502609, rho_c=337.628084504206,
           eps_f=0.1299725503, eps_r=0.1627287372, poro=0.7072987125)


def variant(xx, yy):
    """XX en 0.01 g/cm³, YY en % masse de résine -> bilan de phase solide."""
    m_f = 10.0 * xx                       # kg/m³ : la densité de la préforme
    if not 0.0 < yy < 100.0:
        raise ValueError("YY doit être dans ]0, 100[")
    m_r = m_f * yy / (100.0 - yy)         # kg/m³ de résine vierge
    m_rc = m_r * CHAR_YIELD               # kg/m³ de résine charbonnée

    eps_f = m_f / RHO_FIBER
    eps_r = m_r / RHO_RESIN
    eps_rc = m_rc / RHO_RESIN_CHAR

    rho_v, rho_c = m_f + m_r, m_f + m_rc
    return dict(
        xx=xx, yy=yy, m_f=m_f, m_r=m_r, m_rc=m_rc,
        rho_v=rho_v, rho_c=rho_c,
        eps_f=eps_f, eps_r=eps_r, eps_rc=eps_rc,
        poro_v=1.0 - eps_f - eps_r,
        poro_c=1.0 - eps_f - eps_rc,
        k=(rho_v - rho_c) / rho_c,                   # = B'g / B'c
        gas=rho_v - rho_c,
        tga_loss=(yy / 100.0) * (1.0 - CHAR_YIELD),  # perte ATG du COMPOSITE
        w_resin=yy / 100.0,
    )


def yy_max(xx):
    """Teneur en résine maximale : tous les pores de la préforme remplis."""
    eps_f = 10.0 * xx / RHO_FIBER
    m_r = (1.0 - eps_f) * RHO_RESIN
    return 100.0 * m_r / (10.0 * xx + m_r)


def parse(arg):
    a, b = arg.replace("-", "/").split("/")
    return float(a), float(b)


# ---------------------------------------------------------------------------

def main():
    args = sys.argv[1:] or ["14/40", "18/50", "18/80"]
    names = [parse(a) for a in args]

    print()
    print("=" * 79)
    print("  VARIANTES ZURAM XX/YY — fibres Calcarb et résine phénolique INCHANGÉES")
    print("=" * 79)
    print(f"  rho_fibre  {RHO_FIBER:9.2f} kg/m³      rendement char "
          f"{CHAR_YIELD:.4f}")
    print(f"  rho_résine {RHO_RESIN:9.2f} kg/m³      (densité de matrice "
          f"supposée inchangée)")

    # --- 1. contrôle de la paramétrisation sur le vrai 18/50 ---------------
    print("\n" + "-" * 79)
    print("  1. CONTRÔLE — le modèle reproduit-il le ZURAM_official ?")
    print("-" * 79)
    nom = variant(18, 50)
    print(f"  {'':22} {'18/50 nominal':>15} {'ZURAM_official':>15} {'écart':>10}")
    for lab, key, ref in (("rho vierge  [kg/m³]", "rho_v", REF["rho_v"]),
                          ("rho char    [kg/m³]", "rho_c", REF["rho_c"]),
                          ("porosité vierge", "poro_v", REF["poro"])):
        v = nom[key]
        print(f"  {lab:22} {v:15.2f} {ref:15.2f} {100*(v-ref)/ref:9.1f} %")
    print(f"\n  L'écart vient de l'anomalie n° 2 : la préforme réelle pèse "
          f"{REF['eps_f']*RHO_FIBER:.0f} kg/m³,")
    print(f"  non les 180 nominaux, et la résine y fait 51.07 % et non 50.")
    print(f"  -> les valeurs ci-dessous sont NOMINALES ; compter ±15 % sur les")
    print(f"     densités absolues tant que la préforme n'est pas pesée.")

    # --- 2. tableau des variantes -----------------------------------------
    print("\n" + "-" * 79)
    print("  2. CE QUE CHAQUE VARIANTE DONNE")
    print("-" * 79)
    rows = [variant(x, y) for x, y in names]
    hdr = [f"{int(r['xx'])}/{int(r['yy'])}" for r in rows]
    def line(lab, key, fmt="{:12.1f}", scale=1.0):
        print(f"  {lab:26}" + "".join(fmt.format(r[key]*scale) for r in rows))
    print(f"  {'':26}" + "".join(f"{h:>12}" for h in hdr))
    print("  " + "-" * 75)
    line("masse fibres    [kg/m³]", "m_f")
    line("masse résine    [kg/m³]", "m_r")
    print()
    line("rho VIERGE      [kg/m³]", "rho_v")
    line("rho CHAR        [kg/m³]", "rho_c")
    print()
    line("porosité vierge", "poro_v", "{:12.4f}")
    line("porosité char", "poro_c", "{:12.4f}")
    line("eps_fibres", "eps_f", "{:12.4f}")
    line("eps_résine", "eps_r", "{:12.4f}")
    print()
    line("gaz libéré      [kg/m³]", "gas")
    line("perte ATG composite [%]", "tga_loss", "{:12.2f}", 100.0)
    line("couplage k = B'g/B'c", "k", "{:12.4f}")
    print()
    print(f"  {'YY max (pores pleins)':26}"
          + "".join(f"{yy_max(r['xx']):11.1f} %" for r in rows))
    print(f"  {'pores remplis de résine':26}"
          + "".join(f"{100*r['eps_r']/(1-r['eps_f']):11.1f} %" for r in rows))
    ck = lambda e: e ** 3 / (1.0 - e) ** 2
    print(f"  {'Carman-Kozeny e^3/(1-e)^2':26}"
          + "".join(f"{ck(r['poro_v']):12.4f}" for r in rows))
    print("\n  La dernière ligne n'est qu'un ORDRE DE GRANDEUR de perméabilité")
    print("  relative (K proportionnel à eps^3/(1-eps)^2), pas une mesure. Le")
    print("  classeur donne 1.06e-11 m2 pour le ZURAM vierge mesuré.")

    # --- 3. séparation des rôles ------------------------------------------
    print("\n" + "-" * 79)
    print("  3. SÉPARATION DES RÔLES  —  XX est une échelle, YY est la physique")
    print("-" * 79)
    print("  Le couplage k ne dépend QUE de YY. En posant w = YY/100 :")
    print()
    print("      k = (rho_v - rho_c)/rho_c = w(1-Y) / (1 - w(1-Y))     Y = rendement char")
    print()
    print("  XX disparaît de l'expression : il se simplifie. Vérification —")
    print(f"  {'':14}" + "".join(f"{'XX=' + str(x):>12}" for x in (10, 14, 18, 25, 40)))
    for yy in (40, 50, 80):
        w = yy / 100.0
        formule = w * (1 - CHAR_YIELD) / (1 - w * (1 - CHAR_YIELD))
        vals = "".join(f"{variant(x, yy)['k']:12.6f}" for x in (10, 14, 18, 25, 40))
        print(f"  YY={yy:<3} k =" + vals + f"   formule {formule:.6f}")
    print("\n  Idem pour la perte de masse en ATG du composite : w(1-Y), sans XX.")
    print("  En revanche la POROSITÉ dépend des deux, et les densités varient")
    print("  proportionnellement à XX à YY fixé.")

    # --- 3bis. ATG ---------------------------------------------------------
    print("\n" + "-" * 79)
    print("  3bis. LES ATG — même courbe, autre amplitude")
    print("-" * 79)
    sf = sum(f for f, _, _, _ in KINETICS)
    print(f"  Cinétique du classeur ('ZURAM_official'!A79:F83), somme des f =")
    print(f"  {sf:.6f}, soit un rendement en char de {1-sf:.4f} pour la RÉSINE.\n")
    print(f"  {'réaction':>9} {'f':>12} {'log10 A':>9} {'E/R [K]':>11} {'m':>6}")
    for i, (f, la, er, m) in enumerate(KINETICS, 1):
        print(f"  {i:>9} {f:12.6f} {la:9.2f} {er:11.1f} {m:6.2f}")

    print("\n  ATG SUR LA RÉSINE SEULE — strictement IDENTIQUE pour toutes les")
    print("  variantes. Mêmes A, E/R, m : mêmes températures de pic, mêmes")
    print("  vitesses, même perte finale de 38.02 %. C'est le même polymère.")
    print("\n  ATG SUR LE COMPOSITE — même FORME, amplitude proportionnelle à YY :")
    print(f"\n  {'variante':>10} {'w_résine':>10} {'perte totale':>14} "
          f"{'masse résiduelle':>18}")
    for r in rows:
        print(f"  {int(r['xx'])}/{int(r['yy']):<7} {r['w_resin']:10.3f} "
              f"{100*r['tga_loss']:13.2f} % {100*(1-r['tga_loss']):17.2f} %")
    print("\n  Les fibres Calcarb ne pyrolysent pas : elles ne font que diluer")
    print("  le signal. Une ATG de composite ne distingue donc PAS deux résines")
    print("  différentes d'une même résine en proportion différente — il faut")
    print("  normaliser par w_résine avant toute comparaison.")

    print("\n  Pour un code qui attend des f sur le COMPOSITE, multiplier :")
    print(f"  {'variante':>10}" + "".join(f"{'f' + str(i):>12}" for i in range(1, 5)))
    for r in rows:
        print(f"  {int(r['xx'])}/{int(r['yy']):<7}"
              + "".join(f"{f * r['w_resin']:12.6f}" for f, _, _, _ in KINETICS))

    print("\n  RÉSERVE PHYSIQUE. Tout ceci suppose une résine identiquement")
    print("  RÉTICULÉE. À 80 % en masse, la résine remplit 62 % du volume")
    print("  poreux : on n'a plus un film mince autour des fibres mais des")
    print("  amas épais, dont la cuisson (exothermie, évacuation de l'eau et")
    print("  de l'ammoniac de la HMTA) n'a aucune raison d'être la même.")
    print("  Sykes montre justement que le degré de cuisson gouverne la")
    print("  composition (composés piégés) et donc le rendement en char.")
    print("  À l'ATG s'ajoutent des effets de transport interne : craquage")
    print("  secondaire des volatils sur un trajet plus long, que")
    print("  Torres-Herrador quantifie par les nombres Bi, PyI et PyII.")
    print("  -> 14/40 est une extrapolation sûre ; 18/80 demande une ATG neuve.")

    # --- 4. la table B' ----------------------------------------------------
    print("\n" + "-" * 79)
    print("  4. LA TABLE B' EST INCHANGÉE")
    print("-" * 79)
    print("  `bprime` ne consomme que trois compositions élémentaires plus")
    print("  T, P et B'g (src/apps/bprime.cpp:304-330) :")
    print()
    print(f"     bord de couche limite  air        O:0.21, N:0.79")
    print(f"     gaz de pyrolyse        résine     "
          + ", ".join(f"{e}:{v}" for e, v in PYRO_GAS.items()))
    print(f"     char                   fibres+résine carbonisée   C:1.0")
    print()
    print("  Aucune des trois ne dépend de XX ni de YY :")
    print("    - le gaz est produit par la RÉSINE SEULE, identique dans toutes")
    print("      les variantes ; XX/YY changent la QUANTITÉ produite, pas la")
    print("      composition — et cette quantité est justement B'g, le paramètre")
    print("      d'entrée de la table ;")
    print("    - fibres de carbone et résine carbonisée sont toutes deux du")
    print("      carbone pur, donc C:1.0 quelles que soient les proportions.")
    print()
    print("  Densité, porosité et fraction volumique sont des grandeurs")
    print("  VOLUMIQUES, alors que le bilan de surface est écrit sur des")
    print("  fractions massiques élémentaires, donc normalisées.")
    print()
    print("  -> une seule table B' sert TOUTES les variantes XX/YY.")
    print("     Ce qui change est le POINT DE FONCTIONNEMENT sur cette table :")
    print("     la droite B'g = k·B'c, dont la pente k est donnée ci-dessus.")

    # --- 5. vitesse de recul ----------------------------------------------
    print("\n" + "-" * 79)
    print("  5. ATTENTION — même B'c ne veut pas dire même récession")
    print("-" * 79)
    print("  s_dot = m_c/rho_c = B'c · m_e / rho_c : la vitesse de recul est")
    print("  inversement proportionnelle à rho_c.")
    print()
    ref_c = variant(18, 50)["rho_c"]
    print(f"  {'variante':>12} {'rho_char':>12} {'recul relatif à 18/50':>24}")
    for r in rows:
        print(f"  {int(r['xx'])}/{int(r['yy']):<10} {r['rho_c']:12.1f} "
              f"{ref_c/r['rho_c']:20.2f} x")

    print("\n" + "=" * 79)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
