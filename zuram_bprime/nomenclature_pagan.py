#!/usr/bin/env python
"""
La nomenclature « ZURAM XX/YY » confrontée aux variantes réellement fabriquées.

`resine_zuram.md` § 1 posait que XX est la densité de la préforme et que YY est
la teneur en résine en % masse — le premier point établi par le classeur, le
second seulement DÉDUIT, faute de source. La table 1 de

  [P] A. S. Pagan, B. Massuti-Ballester, G. Herdrich, "The Ablation
      Performance and Dynamics of the Heat Shield Material ZURAM®",
      Institut für Raumfahrtsysteme, Universität Stuttgart, 2017.

tranche la question : elle donne cinq variantes avec leurs densités, dont deux
obtenues en réduisant délibérément la teneur en résine.

Trois vérifications :

  1. LE NOM. « resin content reduced by 25 % » -> 18/43 et « by 50 % » -> 18/33.
     Si YY est la fraction massique de résine, l'arithmétique doit tomber juste.
  2. LES DENSITÉS MESURÉES. Écart au nominal, et cas des variantes dont la
     préforme change (18/50-3 au SiC, M/50 à la mullite).
  3. LE RENDEMENT EN CHAR. Celui qu'impliquent les densités de [P], comparé au
     0.6198 du classeur VKI.

Usage :
    python nomenclature_pagan.py
"""

from zuram_variantes import CHAR_YIELD, parse, variant

# [P] table 1 — variante, rho_pref, rho_comp, rho_char [kg/m³], commentaire
PAGAN = [
    ("18/50",   180.0, 370.0, 265.0, "standard quality"),
    ("18/50-3", 350.0, 470.0, 435.0, "préforme carbone + revêtement CVI-SiC"),
    ("18/43",   180.0, 340.0, 250.0, "resin content reduced by 25%"),
    ("18/33",   180.0, 280.0, 225.0, "resin content reduced by 50%"),
    ("M/50",    128.0, 360.0, 230.0, "préforme mullite au lieu de carbone"),
]
PREFORM_18 = 180.0     # kg/m³, Calcarb CBCF 18/2000 ('Calcarb_official'!D15)


def yy_of(name):
    """42.857 -> 43 : la fraction massique de résine, arrondie à l'entier."""
    return name.split("/")[1].split("-")[0]


# ---------------------------------------------------------------------------

def check_naming():
    print("=" * 78)
    print("1. LE NOM SUIT-IL LA TENEUR EN RÉSINE ?")
    print("=" * 78)
    print("  Hypothèse : préforme 180 kg/m³ et, pour 18/50, 50 % de résine en")
    print("  masse — donc 180 kg/m³ de résine et 360 kg/m³ de composite.")
    print("  [P] décrit 18/43 et 18/33 comme la MÊME recette à teneur réduite.\n")
    print(f"  {'réduction':>10} {'résine':>9} {'composite':>11} "
          f"{'w = résine/composite':>21} {'nom prédit':>12} {'nom réel':>10}")
    ok = True
    for red, real in ((0.00, "18/50"), (0.25, "18/43"), (0.50, "18/33")):
        r = PREFORM_18 * (1.0 - red)
        comp = PREFORM_18 + r
        w = 100.0 * r / comp
        pred = f"18/{round(w):02d}"
        ok &= (pred == real)
        print(f"  {-100*red:9.0f}% {r:9.1f} {comp:11.1f} {w:20.3f} % "
              f"{pred:>12} {real:>10}")
    print(f"\n  -> {'LES TROIS NOMS TOMBENT EXACTEMENT' if ok else 'ÉCHEC'}")
    print("  L'arithmétique est sans ambiguïté : 180·0.75/(180+135) = 42.857 %")
    print("  et 180·0.50/(180+90) = 33.333 %. Aucune autre lecture de YY ne")
    print("  produit à la fois 50, 43 et 33 à partir de -0 %, -25 % et -50 %.")
    return ok


def check_measured():
    print("\n" + "=" * 78)
    print("2. LES DENSITÉS MESURÉES")
    print("=" * 78)
    print(f"  {'variante':>9} {'ρ_pref':>8} {'ρ_comp':>8} {'résine':>8} "
          f"{'w mesuré':>10} {'YY du nom':>10} {'écart':>8}")
    for name, rp, rc, _, _ in PAGAN:
        r = rc - rp
        w = 100.0 * r / rc
        yy = float(yy_of(name))
        print(f"  {name:>9} {rp:8.0f} {rc:8.0f} {r:8.0f} {w:9.2f} % "
              f"{yy:10.0f} {w - yy:+7.2f}")
    print("\n  Les trois variantes à préforme carbone standard se tiennent à")
    print("  +1.4, +4.1 et +2.7 points au-dessus de leur nom : la densité")
    print("  réelle dépasse systématiquement le nominal — exactement ce que")
    print("  le DLR explique par l'amélioration de l'infiltration (moule RTM),")
    print("  cf. `resine_zuram.md` § 3.")
    print("\n  Les deux autres NE SUIVENT PAS la règle, et c'est instructif :")
    print("    18/50-3  w = 25.5 % — le revêtement CVI-SiC alourdit la préforme")
    print("             de 180 à 350 kg/m³ ; la résine, elle, est inchangée.")
    print("    M/50     w = 64.4 % — la préforme mullite ne pèse que 128 kg/m³")
    print("             et laisse plus de volume à imprégner.")
    print("  Dans les deux cas le « 50 » est celui de la RECETTE de base, repris")
    print("  tel quel. YY nomme donc la formulation visée, pas le composite tel")
    print("  que construit. Et le premier champ nomme la préforme — d'où le")
    print("  « M » de M/50, qui remplace le nombre par une lettre.")


def check_char_yield():
    print("\n" + "=" * 78)
    print("3. LE RENDEMENT EN CHAR — UN DÉSACCORD À SIGNALER")
    print("=" * 78)
    name, rp, rc, rch, _ = PAGAN[0]
    r_v, r_c = rc - rp, rch - rp
    y = r_c / r_v
    print(f"  [P] mesure ρ_char pour le 18/50 seul, puis met les autres à")
    print(f"  l'échelle « assuming identical ratios of solid resin residue ».\n")
    print(f"  18/50 : résine vierge {r_v:.0f} -> résine charbonnée {r_c:.0f} kg/m³")
    print(f"          rendement en char = {r_c:.0f}/{r_v:.0f} = {y:.4f}\n")
    print("  Contrôle de la mise à l'échelle annoncée par la légende :")
    print(f"  {'variante':>9} {'prédit':>9} {'table [P]':>11} {'écart':>8}")
    for nm, p, c, ch, _ in PAGAN:
        if nm.startswith("18/") and nm != "18/50-3":
            pred = p + (c - p) * y
            print(f"  {nm:>9} {pred:9.1f} {ch:11.0f} {pred - ch:+8.1f}")
    print("  -> la légende est cohérente avec ses propres chiffres.\n")
    print(f"  MAIS le classeur VKI donne {CHAR_YIELD:.4f}, pas {y:.4f}.")
    pred = rp + r_v * CHAR_YIELD
    print(f"  Avec {CHAR_YIELD:.4f} : ρ_char(18/50) = {rp:.0f} + {r_v:.0f}·"
          f"{CHAR_YIELD:.4f} = {pred:.1f} kg/m³,")
    print(f"  alors que [P] mesure {rch:.0f} — soit {pred - rch:.0f} kg/m³ de trop,")
    print(f"  {100*(pred-rch)/rch:.0f} % d'écart.\n")
    print("  Explication la plus probable, à confirmer : les deux grandeurs ne")
    print("  décrivent pas le même état. Le 0.6198 du classeur vient d'une ATG")
    print("  arrêtée vers 800-1000 °C ; le char de [P] est celui d'échantillons")
    print("  passés en jet plasma, donc bien au-delà de 2000 K. Or le modèle")
    print("  cinétique lui-même garde une queue algébrique au-delà de l'ATG")
    print("  (cf. `variantes_zuram.md`) : un char plus chaud est plus léger.")
    print("\n  CONSÉQUENCE PRATIQUE. Pour une réponse matériau en ablation, le")
    print("  0.447 de [P] est sans doute plus représentatif que le 0.6198 du")
    print("  classeur. Le couplage k = (ρ_v - ρ_c)/ρ_c en dépend directement :")
    for lab, yy in (("classeur VKI", CHAR_YIELD), ("[P], 18/50 mesuré", y)):
        k = (rc - (rp + r_v * yy)) / (rp + r_v * yy)
        print(f"      {lab:22s} rendement {yy:.4f} -> k = {k:.4f}")
    print("  La table B', elle, reste inchangée : elle ne consomme pas le")
    print("  rendement en char (cf. `variantes_zuram.md` § 2).")
    return y


def check_invented():
    print("\n" + "=" * 78)
    print("4. NOS VARIANTES INVENTÉES, REPLACÉES DANS LA FAMILLE RÉELLE")
    print("=" * 78)
    fam = sorted(float(yy_of(n)) for n, *_ in PAGAN if n.startswith("18/"))
    print(f"  Famille réellement fabriquée à préforme carbone : YY = "
          f"{', '.join(f'{v:.0f}' for v in fam)}")
    print(f"  -> le DLR n'a exploré que la RÉDUCTION de résine, jamais")
    print(f"     l'augmentation. Plage réelle : {min(fam):.0f} à {max(fam):.0f} %.\n")
    for lab in ("14/40", "18/80"):
        xx, yy = parse(lab)
        v = variant(xx, yy)
        inside = min(fam) <= yy <= max(fam)
        print(f"  {lab} : ρ_vierge {v['rho_v']:.1f} kg/m³, "
              f"porosité {v['poro_v']:.3f}")
        if lab == "14/40":
            print(f"      YY = 40 est DANS la plage explorée (33-50) ->")
            print(f"      extrapolation soutenue par 18/43 et 18/33.")
            print(f"      Mais XX = 14 suppose une préforme Calcarb CBCF 14,")
            print(f"      absente de la famille : toutes gardent CBCF 18.")
        else:
            print(f"      YY = 80 est HORS de la plage explorée (33-50), et")
            print(f"      dans la direction que le DLR n'a jamais prise.")
            print(f"      La réserve de `variantes_zuram.md` § 4 est donc")
            print(f"      confirmée empiriquement, pas seulement théorique.")


def main():
    print()
    print("  NOMENCLATURE ZURAM XX/YY — CONFRONTATION À PAGAN et al. (2017)")
    print()
    ok = check_naming()
    check_measured()
    check_char_yield()
    check_invented()
    print("\n" + "=" * 78)
    print(f"  [{'OK' if ok else 'ÉCHEC'}] YY = teneur en résine en % masse "
          f"de la recette — CONFIRMÉ")
    print("  [!] rendement en char : 0.4474 d'après [P] contre 0.6198 au")
    print("      classeur — désaccord réel, à trancher avant tout calcul de")
    print("      réponse matériau.")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
