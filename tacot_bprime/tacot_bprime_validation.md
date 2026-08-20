# Validation de la table B' TACOT contre la référence TACOT 3.0

## 1. Source de référence

Classeur `TACOT_3.0_1.xls`, feuille **`B' table`**.

> *B' table (0.001, 0.01, 0.1 and 1 atm, air) — initially computed with Mutation
> using the thermophysical database of CEA, here the results are calculated with
> TARGET.*
> *Hypotheses: equal diffusion coefficients, equilibrium chemistry at the wall,
> air: O2/N2 = 0.21/0.79 molar.*

- [1] J. de Muelenaere, J. Lachaud, N. N. Mansour, T. E. Magin, *Stagnation line
  approximation for ablation thermochemistry*, 42nd AIAA Thermophysics
  Conference, 2011.
- [2] D. Bianchi, *A CEA-based Chemical Equilibrium Solver for Gas/Surface
  Thermochemistry and Thermochemical Tables Generation*, CRAS-TTG-1001,
  Sapienza University of Rome, 2013.

La table contient **15 100 points** : 4 pressions × 25 valeurs de B'g × 151
températures.

| Grandeur | Valeurs |
|----------|---------|
| Pression | 0.001, 0.01, 0.1, 1.0 atm (stockées en bar) |
| B'g | 0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0, 4.0, 5.5, 7.5, 10.0 |
| Température | 250 → 4000 K, pas 25 K |
| Colonnes | `p (bar)`, `B'g`, `Temp (K)`, `B'c`, `Hw (kJ/kg)` |

> Attention aux unités : la pression est en **bar** (`1.01325D-03` = 0.001 atm,
> notation Fortran `D`) et l'enthalpie en **kJ/kg**, alors que `bprime` attend
> des Pa et sort des MJ/kg.

---

## 2. Mise en données corrigée

### 2.1 Gaz de pyrolyse

La composition retenue précédemment (phénol pur, `C:6, H:6, O:1`) **ne
correspond pas** au TACOT. Le classeur donne, feuille `Pyrolysis model` :

```
C : 0.206     H : 0.679     O : 0.115      (fractions molaires élémentaires)
```

Elle dérive de la composition moléculaire de la feuille
`Pyrolysis-gas chemistry` (Sykes, *NASA TN D-3810*, 1967) :

| Espèce | H2 | H2O | CH4 | C6H5OH | CO | CO2 | C6H6 |
|--------|-----|-----|-----|--------|-----|-----|------|
| x [-] | 0.4992 | 0.2336 | 0.1000 | 0.0891 | 0.0576 | 0.0157 | 0.0048 |

Vérification : C = 0.7368, H = 2.4291, O = 0.4117 mol → normalisé
**C:0.2059, H:0.6790, O:0.1151**, soit exactement les valeurs tabulées.

Le gaz réel est donc bien plus **riche en hydrogène** que le phénol pur
(H2 représente à lui seul la moitié des moles), ce qui change nettement
h_w et B'c dès que B'g > 0.

### 2.2 Liste d'espèces

La référence utilise **25 espèces gazeuses**, obtenues par réduction d'une liste
CEA de 112 espèces :

```
C H O N CH4 CN CO CO2 C2 C2H C2H2,acetylene C3 C4 C4H2,butadiyne C5
HCN H2 H2O N2 CH2OH CNN CNC CNCOCN C6H6 HNC
```

Le fichier `data/mixtures/tacot-air_25.xml` reproduit cette liste **à
l'identique**, plus `C(gr)` pour la phase condensée du char. `C(gr)` ne fait pas
partie de la liste *gaz* de la référence, mais le solveur multiphase en a besoin
pour que le carbone du char en excès puisse se condenser ; sans lui,
Y_w,C → 1 et B'c sature sur une valeur non physique.

### 2.3 Char

Char = carbone pur (`C:1.0`), hypothèse standard de la formulation B'.

---

## 3. Résultats

Commande :

```bash
python tacot_bprime_validation.py --mixture tacot-air_25
```

### 3.1 Découpage en régimes

Trois zones sont distinguées, car toutes ne sont pas comparables :

| Régime | Définition | Points |
|--------|-----------|--------|
| **ACTIF** | référence non nulle et hors palier → comparable | 7 609 |
| **ZERO** | référence exactement nulle (le char n'est pas consommé) | 5 326 |
| **PALIER** | référence gelée après sublimation complète → non comparable | 2 165 |

Le **palier** est un artefact du générateur de la référence : au-delà de la
sublimation totale il n'existe plus de phase condensée, et TARGET répète la
dernière valeur convergée de B'c *et* de h_w jusqu'à 4000 K. Mutation++ continue
au contraire à résoudre l'équilibre et sature sur sa propre limite (quantité
finie de char ajoutée, `max(200, 100·B'g)`). Les deux valeurs sont des artefacts
distincts dans cette zone : elle est exclue des statistiques et grisée sur la
carte d'écart.

Températures de début de palier (sublimation complète) — la progression avec la
pression suit Clausius-Clapeyron :

| P [atm] | 0.001 | 0.01 | 0.1 | 1 |
|---------|-------|------|-----|---|
| T [K] (B'g = 0) | 3075 | 3325 | 3625 | 3950 |
| T [K] (B'g = 10) | 3025 | 3275 | 3550 | 3900 |

### 3.2 Accord — régime actif

| Grandeur | Moyenne | Médiane | p95 | Max |
|----------|---------|---------|-----|-----|
| **B'c** | 0.433 % | 0.011 % | 2.49 % | 40.5 % |
| **h_w** | 0.137 % | 0.083 % | 0.41 % | 1.37 % |

- **B'c** : 89.2 % des points sous 1 %, 98.3 % sous 5 %, 99.7 % sous 10 %
- **h_w** : 99.7 % des points sous 1 %, **100 % sous 5 %**

### 3.3 Accord — régime B'c = 0

Sur les 5 326 points où la référence impose B'c = 0, Mutation++ donne
également B'c < 10⁻⁶ sur **5 325 points (99.98 %)**. Le seuil d'apparition de
l'ablation est donc reproduit.

### 3.4 Écart par pression (B'c, régime actif)

| P [atm] | moyenne | médiane | p95 |
|---------|---------|---------|-----|
| 0.001 | 0.180 % | 0.006 % | 0.69 % |
| 0.01 | 0.170 % | 0.008 % | 0.83 % |
| 0.1 | 0.443 % | 0.012 % | 2.65 % |
| 1 | 0.798 % | 0.021 % | 4.79 % |

### 3.5 Où sont les écarts résiduels

Seuls **24 points sur 7 609 (0.32 %)** dépassent 10 % d'écart :

- **20 points à T = 3925 K, P = 1 atm** — soit exactement un pas de grille avant
  le palier de la référence (3950 K). C'est le *genou de sublimation*, où B'c est
  quasi vertical : un décalage de quelques kelvins sur le seuil de sublimation
  produit mécaniquement un grand écart relatif. L'écart y est systématique
  (≈ −11 %) et non aléatoire.
- **3 points vers 700–975 K** où B'c ≈ 5·10⁻⁴ : l'écart relatif est grand mais
  l'écart absolu est négligeable.
- 1 point à 2775 K.

Les deux bandes visibles sur la carte d'écart correspondent aux deux genoux de
la courbe (oxydation vers 500–1000 K, sublimation vers 3000–4000 K), là où la
dérivée dB'c/dT est maximale.

---

## 4. Effet de la liste d'espèces

Le même calcul mené avec `tacot-air_35.xml` (24 espèces de la référence + CH,
CH2, CH3, C2H4, C2H6, C3H3 ×2, C6H5O, C6H5OH, O2, OH, mais sans HNC) est
nettement moins bon :

| Mélange | B'c moyen | B'c médian | < 1 % | h_w moyen |
|---------|-----------|------------|-------|-----------|
| **tacot-air_25** (liste de la référence) | **0.433 %** | **0.011 %** | **89.2 %** | **0.137 %** |
| tacot-air_35 (liste étendue) | 2.409 % | 0.434 % | 61.3 % | 0.298 % |

L'écart se creuse avec la pression (3.5 % en moyenne à 1 atm contre 1.0 % à
0.001 atm) : les espèces supplémentaires, favorisées à haute densité, déplacent
l'équilibre à la paroi. **Pour reproduire la table de référence, il faut
utiliser sa liste d'espèces** — c'est `tacot-air_25` qui doit être employé.

---

## 5. Conclusion

Après correction de la composition du gaz de pyrolyse et alignement de la liste
d'espèces, Mutation++ reproduit la table B' de référence du TACOT avec un écart
médian de **0.011 % sur B'c** et **0.083 % sur h_w**, sur les 7 609 points
comparables de la grille. Les écarts résiduels sont localisés aux genoux de
sublimation et aux valeurs de B'c proches de zéro, et sont cohérents avec la
sensibilité intrinsèque de ces zones.

---

## 6. Fichiers produits

| Fichier | Contenu |
|---------|---------|
| `tacot_bprime_vs_ref_tacot-air_25.csv` | comparaison point par point (15 100 lignes, colonne `regime`) |
| `tacot_bprime_vs_ref_tacot-air_25.png` | superposition calcul / référence pour 4 valeurs de B'g |
| `tacot_bprime_err_tacot-air_25.png` | carte d'écart relatif (palier grisé) |
| `tacot_bprime_vs_ref_tacot-air_35.*` | idem pour la liste étendue |

## 7. Reproduire

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release .
cmake --build build --target bprime -j

cd tacot_bprime
python tacot_bprime_validation.py --mixture tacot-air_25 --ref /chemin/TACOT_3.0_1.xls
```

Dépendances : `xlrd` (lecture du `.xls` BIFF), `numpy`, `matplotlib`.
