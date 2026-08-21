#!/usr/bin/env python
"""
Traçabilité de la résine du ZURAM : de la caractérisation VKI au XML.

Rejoue, à partir des seules données publiées, la chaîne qui mène aux trois
nombres inscrits dans `data/mixtures/zuram-*.xml` :

    <composition name="VKIZuramPyroGas">  C:0.171, H:0.722, O:0.107  </composition>

Sources :
    [ODS] Tacot_Zuram_Calcarb_database_v4.3.1.ods, onglet 'ZURAM_official'
          A. Turchi, F. Torres-Herrador, J. Rico Orero (VKI)
    [THo] F. Torres-Herrador, A. Eschenbacher, J. Coheur, J. Blondeau,
          T. E. Magin, K. M. Van Geem, "Decomposition of carbon/phenolic
          composites for aerospace heatshields: Detailed speciation of
          phenolic resin pyrolysis products", Aerospace Science and
          Technology 119 (2021) 107079.

Les six vérifications :

  1. microstructure -> densités vierge et char       ('ZURAM_official'!D12:F21)
  2. rendement en char, par deux voies indépendantes (volumes / cinétique)
  3. enthalpie de formation du vierge                ('ZURAM_official'!D69:D71)
  4. fractions massiques -> molaires -> le XML       ('ZURAM_official'!D87:E89)
  5. fermeture résine - char -> gaz                  (fig. 7 de [THo])
  6. nomenclature « 18/50 »                          ('Calcarb_official'!D15)

Différence de nature avec le TACOT : la composition du gaz du ZURAM est une
donnée ÉLÉMENTAIRE mesurée directement (fractions massiques), là où celle du
TACOT est reconstruite en sommant les atomes d'une spéciation moléculaire.
Le § 4 est donc une simple conversion, exacte ; c'est le § 5 qui apporte la
validation croisée.

Voir `resine_zuram.md` pour le commentaire détaillé.

Usage :
    python resine_zuram_verification.py
"""

# Masses molaires telles qu'employées par la feuille ('ZURAM_official'!F87:F89)
M = {"C": 12.0107, "H": 1.00794, "O": 15.999, "N": 14.007}

# ---------------------------------------------------------------------------
# Données publiées — saisies telles quelles
# ---------------------------------------------------------------------------

# [ODS] 'ZURAM_official'!D12:F13 — fractions volumiques
VF = {"virgin": {"gas": 0.7072987125,
                 "fibers": 0.1299725503,
                 "resin": 0.1627287372},
      "char":   {"gas": 0.769173653368846,
                 "fibers": 0.1299725503,
                 "resin": 0.100853796331154}}

# [ODS] 'ZURAM_official'!D21:F21 — densités intrinsèques [g/cm³]
RHO_FIBER = 1.5772413793
RHO_RESIN_VIRGIN = 1.3150719638
RHO_RESIN_CHAR = 1.3150719638          # identique : perte de VOLUME, pas de densité

# [ODS] 'ZURAM_official'!D17:E17 — densités moyennes [g/cm³]
RHO_V_SHEET = 0.418998084502609
RHO_C_SHEET = 0.337628084504206

# [ODS] 'ZURAM_official'!D69:D71 — enthalpies de formation à 298 K [J/kg]
HF_CHAR = 0.0
HF_RESIN = -2143700.0
HF_VIRGIN_SHEET = -1094878.03634388

# [ODS] 'ZURAM_official'!B80:B83 — fractions de perte de densité de la résine
KINETICS_F = [0.035070093457862, 0.0276869158877858,
              0.0959813084109907, 0.221495327102286]

# [ODS] 'ZURAM_official'!E87:E89 puis D87:D89 — gaz de pyrolyse, source [4]
GAS_MASS_FRACTIONS = {"C": 0.457, "H": 0.162, "O": 0.381}
GAS_MOLE_SHEET = {"C": 0.170941536278466,
                  "H": 0.722071254424386,
                  "O": 0.106987209297148}

# data/mixtures/zuram-air.xml et zuram-pyrogas.xml
XML_COMPOSITION = {"C": 0.171, "H": 0.722, "O": 0.107}

# [THo] fig. 7 — composition élémentaire de la résine vierge, % masse.
# Relevé graphique : barres calibrées sur les graduations de l'axe. La colonne
# micropyrolyse somme à 99.8 %, or l'article précise qu'elle a été renormalisée
# à 100 % — ce qui contrôle la calibration à 0.2 point près.
RESIN_FIG7 = {"elemental_analyzer": {"C": 75.21, "H": 5.75, "N": 1.44, "O": 14.13},
              "micropyrolysis":     {"C": 74.17, "H": 5.27, "N": 1.04, "O": 19.32}}

# [THo] § 4.2.1 — composition du char mesurée à 800 °C, % masse
CHAR_800C = {"C": 96.36, "N": 0.13, "H": 0.08, "O": 3.43}

# [ODS] 'Calcarb_official' — la préforme seule, qui donne le « 18 » du nom
CALCARB_VF_FIBERS = 0.114140773620799     # 'Calcarb_official'!E12
CALCARB_RHO_FIBER = 1577.0                # 'Calcarb_official'!D19
CALCARB_RHO_MEAN = 180.0                  # 'Calcarb_official'!D15 [kg/m³]

# [ODS] 'ZURAM_official'!B76 — fraction massique de résine de l'échantillon TGA
RESIN_MASS_FRACTION_TGA = 0.5417721519

# [ODS] 'Versions_and_issues'!G25 — masses volumiques citées dans le courriel DLR
RHO_DELIVERED = 430.0     # « the final density of the provided virgin ZURAM (18/50) »
RHO_V12 = 366.0           # « the density indicated (ZURAM V12) »


def normalize(d, total=1.0):
    s = sum(d.values())
    return {k: v / s * total for k, v in d.items()}


def mass_to_mole(y):
    n = {e: y[e] / M[e] for e in y}
    return normalize(n)


# ---------------------------------------------------------------------------
# 1. Microstructure -> densités
# ---------------------------------------------------------------------------

def check_microstructure():
    print("=" * 78)
    print("1. MICROSTRUCTURE  ->  DENSITÉS  ('ZURAM_official'!D12:F21)")
    print("=" * 78)
    print(f"  {'':8} | {'gaz':>10} {'fibres':>10} {'résine':>10} | {'somme':>7}")
    for state in ("virgin", "char"):
        v = VF[state]
        print(f"  {state:8} | {v['gas']:10.6f} {v['fibers']:10.6f} "
              f"{v['resin']:10.6f} | {sum(v.values()):7.4f}")

    rv = VF["virgin"]["fibers"] * RHO_FIBER + VF["virgin"]["resin"] * RHO_RESIN_VIRGIN
    rc = VF["char"]["fibers"] * RHO_FIBER + VF["char"]["resin"] * RHO_RESIN_CHAR
    print(f"\n  rho_v = {VF['virgin']['fibers']:.7f}*{RHO_FIBER} "
          f"+ {VF['virgin']['resin']:.7f}*{RHO_RESIN_VIRGIN}")
    print(f"        = {rv:.15f}   feuille {RHO_V_SHEET:.15f}   "
          f"écart {abs(rv - RHO_V_SHEET):.1e}")
    print(f"  rho_c = {rc:.15f}   feuille {RHO_C_SHEET:.15f}   "
          f"écart {abs(rc - RHO_C_SHEET):.1e}")
    print(f"\n  -> l'unité annoncée [kg/m³] est en fait le g/cm³ : "
          f"{1000*rv:.1f} et {1000*rc:.1f} kg/m³")
    print(f"  -> couplage matériau k = (rho_v - rho_c)/rho_c = {(rv-rc)/rc:.4f}")
    print(f"     gaz libéré {1000*(rv-rc):.1f} kg/m³, "
          f"soit {100*(rv-rc)/rv:.1f} % de la masse vierge")
    ok = abs(rv - RHO_V_SHEET) < 1e-12 and abs(rc - RHO_C_SHEET) < 1e-12
    print(f"  -> {'REPRODUCTION EXACTE' if ok else 'ÉCHEC'}")
    return rv, rc, ok


# ---------------------------------------------------------------------------
# 2. Rendement en char, deux voies
# ---------------------------------------------------------------------------

def check_char_yield():
    print("\n" + "=" * 78)
    print("2. RENDEMENT EN CHAR DE LA RÉSINE — deux voies indépendantes")
    print("=" * 78)
    by_volume = (VF["char"]["resin"] * RHO_RESIN_CHAR
                 / (VF["virgin"]["resin"] * RHO_RESIN_VIRGIN))
    by_kinetics = 1.0 - sum(KINETICS_F)
    print(f"  (a) fractions volumiques, densités intrinsèques identiques :")
    print(f"        {VF['char']['resin']:.9f} / {VF['virgin']['resin']:.9f} "
          f"= {by_volume:.6f}")
    print(f"  (b) 1 - somme des f de la cinétique de pyrolyse :")
    print(f"        1 - {sum(KINETICS_F):.6f} = {by_kinetics:.6f}")
    print(f"\n  écart entre les deux voies : {abs(by_volume - by_kinetics):.1e}")
    ok = abs(by_volume - by_kinetics) < 1e-9
    print(f"  -> rendement en char = {by_volume:.4f}   "
          f"{'CONCORDANT' if ok else 'ÉCHEC'}")
    print("\n  Comme dans le TACOT, aucune cellule ne porte le rendement en char :")
    print("  il est implicite, mais les deux lectures coïncident. Noter qu'ici le")
    print("  classeur suppose une densité intrinsèque de matrice INCHANGÉE, donc")
    print("  une pure perte de volume (eps_m 0.1627 -> 0.1009).")
    return by_volume, ok


# ---------------------------------------------------------------------------
# 3. Enthalpie de formation du vierge
# ---------------------------------------------------------------------------

def check_enthalpy(rho_v):
    print("\n" + "=" * 78)
    print("3. ENTHALPIE DE FORMATION DU VIERGE  ('ZURAM_official'!D69:D71)")
    print("=" * 78)
    w_resin = VF["virgin"]["resin"] * RHO_RESIN_VIRGIN / rho_v
    w_fiber = VF["virgin"]["fibers"] * RHO_FIBER / rho_v
    h = w_resin * HF_RESIN + w_fiber * HF_CHAR
    print(f"  fractions massiques du vierge : résine {w_resin:.6f}, "
          f"fibres {w_fiber:.6f}")
    print(f"  h_f = {w_resin:.6f} * ({HF_RESIN:.0f}) + {w_fiber:.6f} * ({HF_CHAR:.0f})")
    print(f"      = {h:.8f} J/kg")
    print(f"  feuille D71 = {HF_VIRGIN_SHEET:.8f} J/kg   écart {abs(h-HF_VIRGIN_SHEET):.1e}")
    ok = abs(h - HF_VIRGIN_SHEET) < 1e-6
    print(f"  -> {'REPRODUCTION EXACTE' if ok else 'ÉCHEC'}  "
          f"(les fibres de carbone ont h_f = 0 par convention)")
    return ok


# ---------------------------------------------------------------------------
# 4. Gaz de pyrolyse : massique -> molaire -> XML
# ---------------------------------------------------------------------------

def check_gas():
    print("\n" + "=" * 78)
    print("4. GAZ DE PYROLYSE  ('ZURAM_official'!E87:E89 -> D87:D89 -> le XML)")
    print("=" * 78)
    print("  Contrairement au TACOT, la donnée primaire n'est PAS une spéciation")
    print("  moléculaire mais directement la composition ÉLÉMENTAIRE MASSIQUE,")
    print("  mesurée par le VKI (source [4] du classeur). Le passage aux")
    print("  fractions molaires se fait avec les masses molaires de la colonne F.\n")
    x = mass_to_mole(GAS_MASS_FRACTIONS)
    print(f"  {'':3} | {'y (masse)':>10} | {'x calculé':>20} | {'x feuille':>20} | {'écart':>8}")
    worst = 0.0
    for e in ("C", "H", "O"):
        d = abs(x[e] - GAS_MOLE_SHEET[e])
        worst = max(worst, d)
        print(f"  {e:>3} | {GAS_MASS_FRACTIONS[e]:10.3f} | {x[e]:20.16f} | "
              f"{GAS_MOLE_SHEET[e]:20.16f} | {d:8.1e}")
    ok1 = worst < 1e-14
    print(f"\n  écart maximal : {worst:.1e}  -> "
          f"{'CONVERSION EXACTE' if ok1 else 'ÉCHEC'}")

    print(f"\n  Arrondi à 3 décimales pour le XML :")
    ok2 = True
    for e in ("C", "H", "O"):
        r = round(GAS_MOLE_SHEET[e], 3)
        ok2 &= abs(r - XML_COMPOSITION[e]) < 1e-9
        print(f"     {e} : {GAS_MOLE_SHEET[e]:.15f} -> {r:.3f}   "
              f"XML {XML_COMPOSITION[e]:.3f}   {'OK' if abs(r-XML_COMPOSITION[e])<1e-9 else 'ÉCHEC'}")

    mm = sum(GAS_MOLE_SHEET[e] * M[e] for e in GAS_MOLE_SHEET)
    print(f"\n  masse molaire du gaz ZURAM : {mm:.4f} g/mol")
    print(f"  (TACOT pour comparaison : "
          f"{0.206*M['C'] + 0.679*M['H'] + 0.115*M['O']:.4f} g/mol)")
    print("  Le gaz du ZURAM est plus léger : plus riche en H, plus pauvre en C.")
    return ok1 and ok2


# ---------------------------------------------------------------------------
# 5. Fermeture résine - char -> gaz
# ---------------------------------------------------------------------------

def gas_from_resin(resin_pct, char_yield, char_pct):
    """resin_pct, char_pct : % masse. char_yield : fraction massique."""
    char = {e: char_yield * 100.0 * p / 100.0 for e, p in char_pct.items()}
    gas = {e: resin_pct[e] - char.get(e, 0.0) for e in resin_pct}
    if any(v < -1e-9 for v in gas.values()):
        return None
    n = {e: max(gas[e], 0.0) / M[e] for e in ("C", "H", "O")}
    return normalize(n)


def check_closure(char_yield):
    print("\n" + "=" * 78)
    print("5. FERMETURE  résine - char -> gaz   (validation croisée)")
    print("=" * 78)
    print("  La composition du gaz vient d'une source VKI ([4]) ; celle de la")
    print("  résine et du char vient de l'article de Torres-Herrador. Deux")
    print("  campagnes distinctes : la fermeture est un vrai contrôle.\n")

    resin = normalize(RESIN_FIG7["elemental_analyzer"], 100.0)
    print(f"  résine (fig. 7, analyseur élémentaire), renormalisée sur CHNO :")
    print("     " + "  ".join(f"{e} {resin[e]:5.2f} %" for e in ("C", "H", "N", "O")))
    n = {e: resin[e] / M[e] for e in resin}
    print(f"     formule brute : C{n['C']/n['O']:.2f} H{n['H']/n['O']:.2f} O1 "
          f"N{n['N']/n['O']:.2f}")

    s = 100.0 - resin["N"]
    ref = {e: 100 * resin[e] / s for e in ("C", "H", "O")}
    print(f"\n  hors azote : C {ref['C']:.2f}  H {ref['H']:.2f}  O {ref['O']:.2f}")
    print(f"  {'motif':>28} | {'C %':>6} {'H %':>6} {'O %':>6} | {'écart':>6}")
    best, best_d = None, 1e9
    for label, f in (("phénol            C6H6O", {"C": 6, "H": 6, "O": 1}),
                     ("novolac linéaire  C7H6O", {"C": 7, "H": 6, "O": 1}),
                     ("novolac réticulé  C7.5H6O", {"C": 7.5, "H": 6, "O": 1})):
        m = {e: v * M[e] for e, v in f.items()}
        t = sum(m.values())
        p = {e: 100 * m[e] / t for e in m}
        d = max(abs(p[e] - ref[e]) for e in ref)
        if d < best_d:
            best, best_d = label, d
        print(f"  {label:>28} | {p['C']:6.2f} {p['H']:6.2f} {p['O']:6.2f} | {d:6.2f}")
    print(f"  -> motif le plus proche : {best.split()[0]} ({best.split()[-1]}), "
          f"écart {best_d:.2f} point de %")

    print(f"\n  Fermeture, rendement en char {char_yield:.4f} :")
    print(f"  {'char':>18} | {'C':>7} {'H':>7} {'O':>7} | {'écart max':>9}")
    results = {}
    for label, cc in (("carbone pur", {"C": 100.0}),
                      ("mesuré 800 °C", CHAR_800C)):
        g = gas_from_resin(resin, char_yield, cc)
        d = max(abs(g[e] - GAS_MOLE_SHEET[e]) / GAS_MOLE_SHEET[e] for e in g) * 100
        results[label] = d
        print(f"  {label:>18} | {g['C']:7.4f} {g['H']:7.4f} {g['O']:7.4f} | {d:7.1f} %")
    print(f"  {'cible : le XML':>18} | {GAS_MOLE_SHEET['C']:7.4f} "
          f"{GAS_MOLE_SHEET['H']:7.4f} {GAS_MOLE_SHEET['O']:7.4f} | {'—':>9}")

    d = results["carbone pur"]
    print(f"\n  Accord à {d:.1f} % avec l'hypothèse de char de carbone pur — bon,")
    print("  compte tenu que la composition de la résine est relevée sur un")
    print("  graphique (±0.5 point) et provient d'une autre campagne que le gaz.")
    print(f"\n  Le char réel n'est pas du carbone pur : à 800 °C Torres-Herrador")
    print(f"  mesure C {CHAR_800C['C']} / H {CHAR_800C['H']} / N {CHAR_800C['N']} "
          f"/ O {CHAR_800C['O']} %, et signale")
    print("  explicitement le contraste avec l'hypothèse usuelle de 100 % C.")
    print("  Elle reste licite en ablation : la surface dépasse largement 2000 K.")
    return d < 8.0


# ---------------------------------------------------------------------------
# 6. Nomenclature « ZURAM 18/50 »
# ---------------------------------------------------------------------------

def check_nomenclature(rho_v):
    print("\n" + "=" * 78)
    print("6. NOMENCLATURE  « ZURAM® 18/50 »")
    print("=" * 78)

    d = CALCARB_VF_FIBERS * CALCARB_RHO_FIBER
    print("  « 18 » — la préforme Calcarb® CBCF 18/2000 :")
    print(f"     {CALCARB_VF_FIBERS:.15f} x {CALCARB_RHO_FIBER:.0f} = {d:.4f} kg/m³")
    print(f"     'Calcarb_official'!D15 = {CALCARB_RHO_MEAN:.0f} kg/m³   "
          f"écart {abs(d - CALCARB_RHO_MEAN):.1e}")
    print("     -> 18 = 0.18 g/cm³, masse volumique apparente NOMINALE de la préforme")
    ok_18 = abs(d - CALCARB_RHO_MEAN) < 1e-3

    print("\n  « 50 » — non énoncé dans les sources. Lecture « teneur en résine » :")
    m_f = VF["virgin"]["fibers"] * RHO_FIBER * 1000.0
    m_r = VF["virgin"]["resin"] * RHO_RESIN_VIRGIN * 1000.0
    readings = [
        ("microstructure ZURAM_official", 100 * m_r / (m_f + m_r)),
        ("échantillon TGA (B76)", 100 * RESIN_MASS_FRACTION_TGA),
        ("V12 0.366 g/cm³, préforme nominale 180",
         100 * (RHO_V12 - CALCARB_RHO_MEAN) / RHO_V12),
        ("livré 0.43 g/cm³, préforme mesurée %.0f" % m_f,
         100 * (RHO_DELIVERED - m_f) / RHO_DELIVERED),
        ("livré 0.43 g/cm³, préforme nominale 180",
         100 * (RHO_DELIVERED - CALCARB_RHO_MEAN) / RHO_DELIVERED),
    ]
    for label, v in readings:
        print(f"     {label:42s} {v:6.2f} %")
    core = [v for _, v in readings[:3]]
    print(f"\n     les trois premières lectures tiennent dans "
          f"[{min(core):.1f} ; {max(core):.1f}] %")

    print("\n  Lecture concurrente « masse volumique visée 0.50 g/cm³ » : EXCLUE.")
    print(f"     le courriel DLR ('Versions_and_issues'!G25) donne "
          f"{RHO_DELIVERED/1000:.2f} g/cm³ pour")
    print(f"     le ZURAM (18/50) livré et {RHO_V12/1000:.3f} g/cm³ pour la version V12 —")
    print("     deux densités différentes sous le MÊME nom : ce n'est donc pas une densité.")

    print("\n  Les fibres du ZURAM_official pèsent %.1f kg/m³, non 180 :" % m_f)
    print("     anomalie ouverte n° 2 du classeur — variabilité de la préforme,")
    print("     jamais pesée par le DLR ; la montée en densité vient du passage")
    print("     à une infiltration en moule RTM.")

    print("\n  -> 18 : ÉTABLI.   50 : DÉDUCTION cohérente, non énoncée dans les sources.")
    return ok_18


# ---------------------------------------------------------------------------

def main():
    print()
    print("  TRAÇABILITÉ DE LA RÉSINE DU ZURAM® 18/50")
    print("  Caractérisation VKI/AblaNTIS  ->  Tacot_Zuram_Calcarb_database")
    print("  ->  data/mixtures/zuram-air.xml, zuram-pyrogas.xml")
    print()
    print("  ZURAM® : ablateur carbone/phénolique produit par le DLR — préforme")
    print("  Mersen Calcarb CBCF 18/2000 imprégnée de résine phénolique,")
    print("  catalysée à l'hexamine.  ('ZURAM_official'!A4 ; Torres-Herrador § 4.2)")
    print()

    rho_v, _, ok1 = check_microstructure()
    char_yield, ok2 = check_char_yield()
    ok3 = check_enthalpy(rho_v)
    ok4 = check_gas()
    ok5 = check_closure(char_yield)
    ok6 = check_nomenclature(rho_v)

    print("\n" + "=" * 78)
    results = [("microstructure -> densités", ok1),
               ("rendement en char, deux voies concordantes", ok2),
               ("enthalpie de formation du vierge", ok3),
               ("massique -> molaire -> le XML", ok4),
               ("fermeture résine/char/gaz", ok5),
               ("nomenclature : le « 18 » retrouvé", ok6)]
    for label, ok in results:
        print(f"  [{'OK' if ok else 'ÉCHEC'}] {label}")
    print("=" * 78)
    return 0 if all(ok for _, ok in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
