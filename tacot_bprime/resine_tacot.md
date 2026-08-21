# La résine du TACOT : identité, composition, traçabilité

Document de traçabilité. Chaque valeur est rattachée soit à une cellule du
classeur `TACOT_3.0.xls`, soit à un tableau de la référence primaire :

> G. F. Sykes, Jr., *Decomposition characteristics of a char-forming phenolic
> polymer used for ablative composites*, NASA TN D-3810, février 1967.

---

## 1. Identité de la résine

| | Valeur | Source |
|---|---|---|
| Nature | phénol-formaldéhyde, **novolac** (stade « B ») | Sykes, § MATERIAL, p. 2 |
| Référence commerciale | **Union Carbide BRP 5549** | Sykes, § MATERIAL, p. 2 |
| Durcisseur | hexaméthylènetétramine (HMTA) prémélangée | Sykes, § MATERIAL, p. 2 |
| Cuisson | 165 °C 1 h sous vide, puis post-cuisson par paliers jusqu'à 149 °C | Sykes, p. 3 |
| Densité vierge | 1200 kg/m³ | `TACOT_3.0.xls`, `Description!A12` |

Cellule `Description!A12` : *« Matrix: ex-novolac/formaldehyde polymer cured at
300 K, virgin density 1200 kg/m³ »*.

Le novolac est formé **en défaut de formaldéhyde** ; la réticulation est
complétée à la cuisson par le formaldéhyde issu de la décomposition de la HMTA
(l'ammoniac libéré sert de catalyseur). D'où la structure de la p. 3 : cycles
phénoliques pontés par des méthylènes –CH₂–.

**Pourquoi cette résine.** Ce n'est pas un choix de modélisation, c'est le
matériau réellement caractérisé par Sykes. Sykes insiste (Introduction, p. 2)
sur le fait que les études antérieures utilisaient des résines *de laboratoire*
et non celles effectivement mises en œuvre dans les boucliers ; il a donc
caractérisé la résine de qualité industrielle des modèles arc-jet, avec le
**même cycle de cuisson**. Le TACOT hérite de ce choix : c'est la seule
phénolique dont on ait à la fois la composition du gaz de pyrolyse ET l'analyse
élémentaire du résidu, mesurées sur le même échantillon.

---

## 2. Composition élémentaire mesurée (Sykes, table II)

Analyse élémentaire du polymère avant chauffage et du résidu, par paliers :

| T (°C) | C | H | N | O | fraction massique résiduelle |
|---|---|---|---|---|---|
| ambiant | **75.60** | **6.12** | **2.35** | **15.93** | **1.00** |
| 400 | 78.45 | 5.42 | 0.42 | 15.71 | 0.930 |
| 600 | 88.00 | 3.59 | — | 8.41 | 0.73 |
| 850 | 92.65 | 0.90 | — | 6.45 | **0.54** |

Remarques de Sykes :

- L'analyse ne dose que C, H et N ; **l'oxygène est le complément à 100 %**
  (teneur en cendres mesurée à 0.04 %, donc négligeable).
- L'azote (2.35 %) provient de la HMTA. Il **n'est pas** dans le squelette du
  polymère : sa fraction résiduelle chute dès les premiers paliers (0.42 % à
  400 °C), signature d'un composé **piégé**, pas d'un azote de structure.
- De l'eau est également piégée lors du moulage (seule espèce évoluée en
  dessous de 300 °C) ; une post-cuisson prolongée à 225 °C l'élimine.

**Formule brute mesurée**, normalisée sur l'oxygène :

```
C6.32 H6.10 O1  (N0.17)
```

Comparaison avec les motifs idéalisés du réseau novolac :

| Motif | C % | H % | O % |
|---|---|---|---|
| **mesuré, renormalisé hors azote** | **77.42** | **6.27** | **16.31** |
| phénol C₆H₆O | 76.57 | 6.43 | 17.00 |
| novolac linéaire C₇H₆O (1 pont CH₂/cycle) | 79.23 | 5.70 | 15.08 |
| novolac réticulé C₇.₅H₆O (1.5 pont CH₂/cycle) | 80.34 | 5.39 | 14.27 |

Le polymère mesuré est **moins carboné que le réseau idéal** et se situe très
près du **phénol C₆H₆O**. Sykes attribue explicitement cet écart aux composés
piégés (eau, dérivés de la HMTA), p. 13. C'est aussi ce que suggérait la
fermeture inverse de `pyrolysis_gas_from_resin.py` — désormais **confirmée par
la mesure directe** et non plus seulement inférée.

---

## 3. Le gaz de pyrolyse : reconstruction exacte

La table I de Sykes donne la distribution des produits de décomposition
(chromatographie en phase gazeuse, paliers de 50 °C entre 100 et 1000 °C), en
pourcentage des moles totales observées.

Sommes des colonnes (calculées sur les entrées individuelles) :

| | CO₂ | CO | C₆H₆ | C₇H₈ | C₆H₅OH | (CH₃)₂C₆H₃OH | CH₄ | H₂O | H₂ | Σ |
|---|---|---|---|---|---|---|---|---|---|---|
| somme réelle | 1.57 | 5.78 | 0.17 | 0.31 | 7.14 | 1.80 | 10.03 | 23.43 | 50.07 | **100.30** |
| total imprimé | 1.6 | 5.5 | 0.2 | 0.3 | 7.1 | 1.8 | 10.0 | 23.4 | 50.1 | 100.0 |

Le classeur applique ensuite la réduction annoncée en
`'Pyrolysis-gas chemistry'!A8` :

> *« Toluene [C7H8] and Xylenol [C8H10O], present in Sykes but not modeled in
> the finite-rate mechanism of April have been transferred to C6H6 and C6H6O,
> respectively, for consistency of the model. »*

soit C₇H₈ → C₆H₆ (0.17 + 0.31 = 0.48) et 2,4-xylénol → C₆H₅OH (7.14 + 1.80 =
8.94), puis renormalise **par 100.30** (et non par les totaux imprimés).

Résultat — reproduction **exacte** (écart < 1e-15) des mole fractions du
classeur `'Pyrolysis-gas chemistry'!B5:H5` :

| espèce | reconstruit | classeur TACOT 3.0 |
|---|---|---|
| CO₂ | 0.0156530408773679 | 0.0156530408773679 |
| CO | 0.0576271186440678 | 0.0576271186440678 |
| C₆H₆ | 0.0047856430707876 | 0.0047856430707876 |
| C₆H₅OH | 0.0891326021934197 | 0.0891326021934197 |
| CH₄ | 0.1000000000000000 | 0.1000000000000000 |
| H₂O | 0.2335992023928215 | 0.2335992023928220 |
| H₂ | 0.4992023928215354 | 0.4992023928215351 |

Sommation des atomes → **C:0.2059, H:0.6790, O:0.1151**, arrondi en
`'Pyrolysis model'!B4:B6` en **C:0.206, H:0.679, O:0.115** — la valeur portée
par les XML `data/mixtures/tacot-*.xml`.

> Noter que le gaz du modèle TACOT est **sans azote**, alors que la résine en
> contient 2.35 %. C'est cohérent : cet azote est un composé piégé, pas un
> constituant, et la table I de Sykes ne recense aucune espèce azotée.

---

## 4. Le rendement en char de 50 %

**Il n'existe aucune cellule « char yield » dans le classeur.** La valeur de
50 % y est implicite, et apparaît deux fois de façon cohérente :

1. **Bilan de microstructure**, `'Pyrolysis model'!A32:A34` :
   ```
   rho_v = eps_f*rho_f + eps_m*rho_m   = 0.1*1600 + 0.1*1200 = 280 kg/m³
   rho_c = eps_f*rho_f + eps_m_c*rho_m = 0.1*1600 + 0.05*1200 = 220 kg/m³
   hyp: intrinsic densities of virgin and charred matrix are equal
        (only volume loss due to shrinkage during the pyrolysis)
   ```
   La matrice passe de `eps_m = 0.10` à `0.05` à densité intrinsèque constante
   → rendement massique `0.05*1200 / (0.10*1200) = **0.50**`.

2. **Cinétique de Goldstein**, `'Pyrolysis model'!B18:C19` :
   ```
   phase A : rho_v = 300, rho_c = 0
   phase B : rho_v = 900, rho_c = 600
   ```
   → matrice `1200 → 600 kg/m³`, soit **0.50** également.

Les deux lectures sont numériquement équivalentes mais **microstructuralement
contradictoires** (perte de volume à densité constante *vs* perte de densité à
volume constant). Elles donnent les mêmes ρ_v = 280 et ρ_c = 220 kg/m³, seules
grandeurs qui comptent en pratique.

### 50 % ou 54 % ?

Sykes **mesure** 0.54 à 850 °C (table II, DTA). Test de fermeture indépendant —
partir de la résine mesurée (table II, ambiant), retrancher un char de carbone
pur, comparer au gaz de la table I :

| rendement | origine | gaz obtenu | écart max / table I |
|---|---|---|---|
| 0.500 | TACOT (microstructure, Goldstein) | C:0.2317 H:0.6600 O:0.1082 | 12.5 % |
| **0.540** | **Sykes table II, mesuré** | **C:0.2028 H:0.6848 O:0.1123** | **2.4 %** |
| 0.534 | ajustement optimal | C:0.2072 H:0.6811 O:0.1117 | 3.0 % |
| — | cible : table I | C:0.2059 H:0.6790 O:0.1151 | — |

**Conclusion.** Le 0.54 mesuré ferme le bilan à 2.4 % près entre deux
techniques expérimentales totalement indépendantes (analyse élémentaire du
résidu *vs* chromatographie des effluents) — validation forte de la cohérence
interne de Sykes. Le 0.50 du TACOT est une **valeur de modèle arrondie**,
choisie pour donner des densités rondes (280 / 220 kg/m³), pas la mesure.

Cet écart est **sans conséquence sur la table B'** : celle-ci ne consomme que
la composition élémentaire du gaz (prise directement de la table I, pas
reconstruite) et celle du char. Le rendement en char n'intervient que dans la
réponse matériau, via le couplage `k = B'g/B'c = (rho_v - rho_c)/rho_c`.

### Le char n'est pas du carbone pur à 850 °C

Sykes mesure 92.65 % C / 0.90 % H / 6.45 % O à 850 °C. L'hypothèse `C:1.0` des
XML reste légitime : en ablation la surface dépasse largement 2000 K, où la
carbonisation est complète, et les fibres ex-cellulose traitées à 2000 K
(`Description!A11`) sont déjà du carbone pur.

---

## 5. Récapitulatif des sources

| Donnée | Valeur | Localisation exacte |
|---|---|---|
| Type de résine | novolac phénol-formaldéhyde + HMTA | Sykes p. 2 ; `Description!A12` |
| Référence | Union Carbide BRP 5549 | Sykes p. 2 |
| Composition résine vierge | C 75.60 / H 6.12 / N 2.35 / O 15.93 (% masse) | Sykes table II |
| Densité résine vierge | 1200 kg/m³ | `Description!A12` |
| Fibres | carbone ex-cellulose, 2000 K, 1600 kg/m³, L 1 mm, ⌀ 10 µm | `Description!A11` |
| Microstructure | eps_f 0.1 / eps_m 0.1 / porosité 0.8 | `Description!A15:A17` |
| Gaz, composition moléculaire | 7 espèces | `'Pyrolysis-gas chemistry'!B5:H5` ← Sykes table I |
| Gaz, composition élémentaire | C:0.206 H:0.679 O:0.115 | `'Pyrolysis model'!B4:B6` |
| Rendement en char (modèle) | 0.50 | implicite : `'Pyrolysis model'!A32:A33` et `B18:C19` |
| Rendement en char (mesuré) | 0.54 à 850 °C | Sykes table II |
| Char, composition à 850 °C | C 92.65 / H 0.90 / O 6.45 | Sykes table II |
| Chaleur de pyrolyse | 293 kJ/kg (350–850 °C) | Sykes, résumé |
