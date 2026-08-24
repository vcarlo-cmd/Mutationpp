#!/usr/bin/env python
"""
Traçabilité de la résine du TACOT : du papier de Sykes au XML de Mutation++.

Rejoue, à partir des seules données publiées, la chaîne qui mène de l'analyse
expérimentale de 1967 aux trois nombres inscrits dans
`data/mixtures/tacot-*.xml` :

    <composition name="tacot_pyro">C:0.206, H:0.679, O:0.115</composition>

Référence primaire :
    G. F. Sykes, Jr., "Decomposition characteristics of a char-forming phenolic
    polymer used for ablative composites", NASA TN D-3810, 1967.
    - table I  : distribution des produits de décomposition (chromatographie)
    - table II : analyse élémentaire du polymère et du résidu (DTA)

Les quatre vérifications :

  1. table I -> mole fractions du classeur    ('Pyrolysis-gas chemistry'!B5:H5)
  2. mole fractions -> composition élémentaire ('Pyrolysis model'!B4:B6)
  3. identification du motif de résine        (table II vs novolac idéalisé)
  4. fermeture croisée résine/char/gaz        (table II vs table I)

Voir `resine_tacot.md` pour le commentaire détaillé.

Usage :
    python resine_tacot_verification.py
"""

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}

# ---------------------------------------------------------------------------
# Données publiées — saisies telles quelles
# ---------------------------------------------------------------------------

# Sykes, table I : colonnes du tableau, une entrée par palier de 50 °C.
# Unité : pourcentage des moles totales observées.
SYKES_TABLE_I = {
    "CO2":     [0.09, 0.11, 0.32, 0.51, 0.26, 0.17, 0.11],
    "CO":      [0.21, 0.44, 0.87, 1.30, 1.19, 0.77, 0.54, 0.26, 0.20],
    "C6H6":    [0.02, 0.06, 0.06, 0.03],
    "C7H8":    [0.08, 0.13, 0.05, 0.04, 0.01],
    "C6H5OH":  [0.46, 0.81, 2.72, 1.62, 0.79, 0.44, 0.21, 0.09],
    "xylenol": [0.13, 0.25, 0.75, 0.38, 0.14, 0.10, 0.05],
    "CH4":     [0.05, 0.15, 0.75, 1.29, 2.61, 2.35, 1.32, 0.83, 0.40, 0.20, 0.08],
    "H2O":     [1.47, 0.75, 0.48, 0.57, 1.28, 3.44, 3.44, 5.42, 3.35, 2.44,
                0.40, 0.26, 0.13],
    "H2":      [0.76, 1.47, 2.18, 3.65, 5.17, 5.88, 6.64, 7.35, 5.88, 4.50,
                3.65, 2.94],
}

# Totaux tels qu'IMPRIMÉS au bas de la table I (arrondis, somme = 100.0)
SYKES_TABLE_I_PRINTED = {
    "CO2": 1.6, "CO": 5.5, "C6H6": 0.2, "C7H8": 0.3, "C6H5OH": 7.1,
    "xylenol": 1.8, "CH4": 10.0, "H2O": 23.4, "H2": 50.1,
}

# Sykes, table II : analyse élémentaire (% masse) et fraction massique résiduelle
SYKES_TABLE_II = [
    #  T(°C)      C      H      N      O    m/m0
    ("ambiant", 75.60, 6.12, 2.35, 15.93, 1.000),
    ("150",     76.08, 5.81, 2.02, 16.09, 0.994),
    ("200",     76.98, 5.82, 1.72, 15.48, 0.979),
    ("300",     77.14, 5.60, 0.78, 16.48, 0.966),
    ("400",     78.45, 5.42, 0.42, 15.71, 0.930),
    ("450",     79.40, 5.36, 0.40, 14.84, 0.902),
    ("500",     81.53, 4.89, 0.00, 13.58, 0.855),
    ("550",     84.22, 4.47, 0.00, 11.31, 0.800),
    ("600",     88.00, 3.59, 0.00,  8.41, 0.730),
    ("650",     90.81, 2.85, 0.00,  6.35, 0.680),
    ("750",     92.31, 1.54, 0.00,  6.15, 0.590),
    ("850",     92.65, 0.90, 0.00,  6.45, 0.540),
]

# TACOT_3.0.xls, 'Pyrolysis-gas chemistry'!B5:H5
WORKBOOK_MOLE_FRACTIONS = {
    "CO2":    0.0156530408773679,
    "CO":     0.0576271186440678,
    "C6H6":   0.0047856430707876405,
    "C6H5OH": 0.0891326021934197,
    "CH4":    0.1,
    "H2O":    0.23359920239282203,
    "H2":     0.49920239282153506,
}

# TACOT_3.0.xls, 'Pyrolysis model'!B4:B6  ==  data/mixtures/tacot-*.xml
WORKBOOK_ELEMENTAL = {"C": 0.206, "H": 0.679, "O": 0.115}

ATOMS = {
    "CO2":    {"C": 1, "O": 2},
    "CO":     {"C": 1, "O": 1},
    "C6H6":   {"C": 6, "H": 6},
    "C6H5OH": {"C": 6, "H": 6, "O": 1},
    "CH4":    {"C": 1, "H": 4},
    "H2O":    {"H": 2, "O": 1},
    "H2":     {"H": 2},
}


def normalize(d):
    t = sum(d.values())
    return {k: v / t for k, v in d.items()}


def elemental(mole_fractions):
    """Somme les atomes d'une composition moléculaire -> fractions molaires."""
    n = {"C": 0.0, "H": 0.0, "O": 0.0}
    for sp, x in mole_fractions.items():
        for e, nu in ATOMS[sp].items():
            n[e] += nu * x
    return normalize(n)


def mass_percent(atom_counts):
    m = {e: n * M[e] for e, n in atom_counts.items()}
    t = sum(m.values())
    return {e: 100 * v / t for e, v in m.items()}, t


# ---------------------------------------------------------------------------
# 1. Table I -> mole fractions du classeur
# ---------------------------------------------------------------------------

def check_table_I():
    print("=" * 78)
    print("1. SYKES TABLE I  ->  'Pyrolysis-gas chemistry'!B5:H5")
    print("=" * 78)

    totals = {k: sum(v) for k, v in SYKES_TABLE_I.items()}
    S = sum(totals.values())

    print(f"  {'colonne':>9} | {'somme réelle':>13} | {'total imprimé':>14}")
    for k in SYKES_TABLE_I:
        print(f"  {k:>9} | {totals[k]:13.2f} | "
              f"{SYKES_TABLE_I_PRINTED[k]:14.1f}")
    print(f"  {'SOMME':>9} | {S:13.2f} | "
          f"{sum(SYKES_TABLE_I_PRINTED.values()):14.1f}")
    print("\n  Les totaux imprimés sont arrondis et somment à 100.0 ; les")
    print("  colonnes elles-mêmes somment à 100.30. C'est cette dernière")
    print("  valeur que le classeur utilise pour renormaliser.")

    # Réduction annoncée en 'Pyrolysis-gas chemistry'!A8 :
    #   toluène -> benzène,  2,4-xylénol -> phénol
    merged = {
        "CO2":    totals["CO2"],
        "CO":     totals["CO"],
        "C6H6":   totals["C6H6"] + totals["C7H8"],
        "C6H5OH": totals["C6H5OH"] + totals["xylenol"],
        "CH4":    totals["CH4"],
        "H2O":    totals["H2O"],
        "H2":     totals["H2"],
    }
    x = {k: v / S for k, v in merged.items()}

    print("\n  Après report C7H8 -> C6H6 et xylénol -> C6H5OH, "
          "renormalisation par 100.30 :\n")
    print(f"  {'espèce':>8} | {'reconstruit':>20} | {'classeur':>20} | {'écart':>9}")
    worst = 0.0
    for k in x:
        d = abs(x[k] - WORKBOOK_MOLE_FRACTIONS[k])
        worst = max(worst, d)
        print(f"  {k:>8} | {x[k]:20.16f} | "
              f"{WORKBOOK_MOLE_FRACTIONS[k]:20.16f} | {d:9.1e}")
    ok = worst < 1e-14
    print(f"\n  écart maximal : {worst:.1e}   -> "
          f"{'REPRODUCTION EXACTE' if ok else 'ÉCHEC'}")
    return x, ok


# ---------------------------------------------------------------------------
# 2. Mole fractions -> composition élémentaire du XML
# ---------------------------------------------------------------------------

def check_elemental(x):
    print("\n" + "=" * 78)
    print("2. COMPOSITION MOLÉCULAIRE  ->  'Pyrolysis model'!B4:B6  (= le XML)")
    print("=" * 78)
    g = elemental(x)
    print(f"  {'':3} | {'reconstruit':>12} | {'classeur / XML':>15} | {'écart':>8}")
    ok = True
    for e in ("C", "H", "O"):
        d = abs(g[e] - WORKBOOK_ELEMENTAL[e])
        ok &= d < 5e-4
        print(f"  {e:>3} | {g[e]:12.4f} | {WORKBOOK_ELEMENTAL[e]:15.3f} | {d:8.1e}")
    mm = sum(g[e] * M[e] for e in g)
    print(f"\n  masse molaire du gaz : {mm:.3f} g/mol")
    print(f"  -> {'COHÉRENT (arrondi à 3 décimales)' if ok else 'ÉCHEC'}")
    return ok


# ---------------------------------------------------------------------------
# 3. Identification du motif de résine
# ---------------------------------------------------------------------------

def check_resin_structure():
    print("\n" + "=" * 78)
    print("3. IDENTIFICATION DU MOTIF  (Sykes table II, ligne 'ambiant')")
    print("=" * 78)
    _, C, H, N, O, _ = SYKES_TABLE_II[0]
    meas = {"C": C, "H": H, "N": N, "O": O}

    n = {e: meas[e] / M[e] for e in meas}
    print(f"  mesuré (% masse)      : C {C:.2f}  H {H:.2f}  N {N:.2f}  O {O:.2f}")
    print(f"  formule brute (sur O) : C{n['C']/n['O']:.2f} H{n['H']/n['O']:.2f}"
          f" O1 N{n['N']/n['O']:.2f}")

    s = 100.0 - N
    ref = {"C": 100 * C / s, "H": 100 * H / s, "O": 100 * O / s}
    print(f"\n  L'azote vient de l'hexaméthylènetétramine et n'est pas dans le")
    print(f"  squelette (il disparaît dès 400 °C). Hors azote :")
    print(f"    C {ref['C']:.2f}  H {ref['H']:.2f}  O {ref['O']:.2f}\n")

    print(f"  {'motif':>44} | {'C %':>6} {'H %':>6} {'O %':>6} | {'écart':>6}")
    print(f"  {'mesuré, hors azote':>44} | {ref['C']:6.2f} {ref['H']:6.2f} "
          f"{ref['O']:6.2f} | {'—':>6}")
    best, best_d = None, 1e9
    for label, f in (("phénol            C6H6O", {"C": 6, "H": 6, "O": 1}),
                     ("novolac linéaire  C7H6O", {"C": 7, "H": 6, "O": 1}),
                     ("novolac réticulé  C7.5H6O", {"C": 7.5, "H": 6, "O": 1})):
        p, _ = mass_percent(f)
        d = max(abs(p[e] - ref[e]) for e in ref)
        if d < best_d:
            best, best_d = label, d
        print(f"  {label:>44} | {p['C']:6.2f} {p['H']:6.2f} {p['O']:6.2f} "
              f"| {d:6.2f}")
    print(f"\n  -> motif le plus proche : {best.split()[0]} "
          f"({best.split()[-1]}), écart {best_d:.2f} point de %")
    print("  Le polymère mesuré est MOINS carboné que le réseau idéal :")
    print("  Sykes l'attribue aux composés piégés (eau, dérivés de la HMTA).")
    return best_d < 1.5


# ---------------------------------------------------------------------------
# 4. Fermeture croisée : table II  vs  table I
# ---------------------------------------------------------------------------

def gas_from_resin(char_yield):
    """Résine mesurée - char de carbone pur -> gaz (fractions molaires C/H/O)."""
    _, C, H, N, O, _ = SYKES_TABLE_II[0]
    virgin = {"C": C, "H": H, "O": O}          # base 100 g, azote écarté
    gas_mass = dict(virgin)
    gas_mass["C"] -= 100.0 * char_yield
    if gas_mass["C"] < 0:
        return None
    return normalize({e: gas_mass[e] / M[e] for e in gas_mass})


def check_closure(target):
    print("\n" + "=" * 78)
    print("4. FERMETURE CROISÉE  résine (table II) - char pur  ->  gaz (table I)")
    print("=" * 78)
    print("  Deux techniques expérimentales indépendantes doivent se recouper :")
    print("    - table II : analyse élémentaire du résidu (DTA)")
    print("    - table I  : chromatographie des effluents\n")

    _, _, _, _, _, _ = SYKES_TABLE_II[0]
    measured_yield = SYKES_TABLE_II[-1][5]      # 0.54 à 850 °C

    grid = [0.30 + 0.0001 * i for i in range(4000)]
    def err(cy):
        g = gas_from_resin(cy)
        return 1e9 if g is None else sum((g[e] - target[e]) ** 2 for e in target)
    best = min(grid, key=err)

    cases = [
        (0.50, "TACOT : microstructure + Goldstein"),
        (measured_yield, "Sykes table II, mesuré à 850 °C"),
        (best, "ajustement optimal"),
    ]
    print(f"  {'rendement':>10} | {'origine':<36} | {'C':>7} {'H':>7} {'O':>7}"
          f" | {'écart':>7}")
    for cy, label in cases:
        g = gas_from_resin(cy)
        d = max(abs(g[e] - target[e]) / target[e] for e in target) * 100
        print(f"  {cy:10.3f} | {label:<36} | {g['C']:7.4f} {g['H']:7.4f} "
              f"{g['O']:7.4f} | {d:6.1f} %")
    print(f"  {'—':>10} | {'cible : gaz de la table I':<36} | "
          f"{target['C']:7.4f} {target['H']:7.4f} {target['O']:7.4f} | {'—':>7}")

    g = gas_from_resin(measured_yield)
    d = max(abs(g[e] - target[e]) / target[e] for e in target) * 100
    print(f"\n  Le rendement MESURÉ (0.54) ferme le bilan à {d:.1f} % près entre")
    print("  deux mesures indépendantes : cohérence interne de Sykes validée.")
    print("  Le 0.50 du TACOT est une valeur de MODÈLE, arrondie pour donner")
    print("  des densités rondes :")
    print(f"    rho_v = 0.1*1600 + 0.1*1200      = "
          f"{0.1*1600 + 0.1*1200:.0f} kg/m³")
    print(f"    rho_c = 0.1*1600 + 0.05*1200     = "
          f"{0.1*1600 + 0.05*1200:.0f} kg/m³")
    print("  (Goldstein donne la même chose autrement : matrice 1200 -> 600,")
    print("   soit (0+600)/(300+900) = 0.50.)")
    print("\n  Sans conséquence sur la table B' : celle-ci prend la composition")
    print("  élémentaire du gaz DIRECTEMENT de la table I, sans passer par le")
    print("  rendement en char. Ce dernier n'intervient que dans la réponse")
    print("  matériau, via k = B'g/B'c = (rho_v - rho_c)/rho_c.")
    return d < 5.0


# ---------------------------------------------------------------------------

def main():
    print()
    print("  TRAÇABILITÉ DE LA RÉSINE DU TACOT")
    print("  Sykes, NASA TN D-3810 (1967)  ->  TACOT_3.0.xls  ->  "
          "data/mixtures/tacot-*.xml")
    print()
    print("  Résine : phénol-formaldéhyde novolac Union Carbide BRP 5549,")
    print("  prémélangée à l'hexaméthylènetétramine, cuite à 165 °C puis")
    print("  post-cuite jusqu'à 149 °C.  (Sykes § MATERIAL ; Description!A12)")
    print()

    x, ok1 = check_table_I()
    ok2 = check_elemental(x)
    ok3 = check_resin_structure()
    ok4 = check_closure(elemental(x))

    print("\n" + "=" * 78)
    results = [("table I -> mole fractions du classeur", ok1),
               ("mole fractions -> élémentaire du XML", ok2),
               ("motif de résine identifié", ok3),
               ("fermeture croisée table II / table I", ok4)]
    for label, ok in results:
        print(f"  [{'OK' if ok else 'ÉCHEC'}] {label}")
    print("=" * 78)
    return 0 if all(ok for _, ok in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
