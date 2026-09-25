# Table B' du carbure de silicium (SiC)

**SiC** massif (β-SiC dense, fritté ou CVD) dans l'air, sans pyrolyse
(B'g = 0). Même grille et mêmes fichiers de sortie que le carbone et la
silice : 25 isobares de 10⁻³ à 10³ atm, Tw de 300 à 5000 K par pas de 25 K.

| | |
|--|--|
| mélange | `data/mixtures/sic-air.xml` : 30 espèces, dont 9 phases condensées |
| char | `sic` : Si:1.0, C:1.0 (y_Si = 0.7005, y_C = 0.2995) |
| binaire | `bprime` généralisé : `-char sic -char-elem N` |
| couplage | k = B'g/B'c = 0 (pas de gaz de pyrolyse) |

```bash
cd sic_bprime
export MPP_DATA_DIRECTORY=../data     # si ce n'est pas déjà fait
python sic_bprime.py                   # ~1 h (4725 points, calculés un par un)
```

| fichier | contenu |
|--|--|
| `sic_bprime_table.csv` | table complète : P_atm, Tw, B'c, h_w [MJ/kg], fractions molaires des 30 espèces |
| `sic_bprime_bc_table.csv` | `Tw_K, P_bar, Bc` (format long nTw × nP, unités SI) |
| `sic_bprime_hw_table.csv` | `Tw_K, P_bar, hw_Jkg` (format long nTw × nP, unités SI) |
| `sic_bprime_phases.csv` | phases condensées présentes en chaque point, méthode de calcul, contrôle de stabilité |
| `sic_bprime_table.png` | B'c et h_w, isobares puissances de 10 |

---

## 1. Pourquoi `-char-elem N` et pas `Si`

`bprime` écrit le bilan de masse sur un seul élément k :

```
B'c = (Y_e,k − Y_w,k) / (Y_w,k − Y_c,k)
```

Ce bilan suppose que le char passe en bloc dans le gaz. C'est toujours vrai
pour le graphite et quasiment toujours vrai pour la silice, mais **pas pour le SiC**.
Selon (Tw, P), le SiC laisse à la paroi de la silice (oxydation passive), du
silicium liquide ou un résidu de carbone (décomposition
SiC → Si(g) + C(gr)). Un suivi sur Si ou sur C donne alors un B'c faux :

| point | phases à la paroi | suivi N | suivi Si | suivi C |
|--|--|--:|--:|--:|
| 2100 K, 1 atm | SiC | 0.2920 | 0.2920 | 0.2920 |
| 2600 K, 1 atm | SiC + Si(L) | **0.258** | 0.151 | 0.611 |
| 2500 K, 0.001 atm | C(gr) | **185** | **0** | 4.3 |

L'azote de l'air est inerte et n'entre dans aucune phase condensée. Son
bilan donne la masse **nette** gazéifiée par la paroi, celle qu'attend le
bilan d'énergie de surface :

```
B'c = Y_e,N / Y_w,N − 1
```

Pour le carbone, ce suivi redonne exactement `carbon_bprime_table.csv`.
Pour la silice, il redonne `silice_bprime_table.csv` à ~1 % près ; l'écart
vient des masses molaires codées en dur dans `bprime_silica`. En dessous de
~1e-12, B'c(N) n'est que du bruit d'arrondi.

## 2. Pourquoi le calcul est fait point par point

Le solveur multiphase de Mutation++ boucle indéfiniment sur ~7 % des points de
la grille. Tous se situent aux frontières de phases : SiO2/SiO, SiC/C(gr),
SiC/Si(L). Il y a deux causes. Une phase ajoutée à 1e-6 mole est retirée au
pas de Newton suivant puis ré-ajoutée, en boucle. Ailleurs, le pas de
continuation décroît sans fin. `sic_bprime.py` procède donc ainsi :

1. Il appelle `bprime` sur **un seul point** (Tw, P), avec un délai de 2 s.
2. Il vérifie le résultat par un critère indépendant (polynômes NASA-9 lus
   dans `nasa9.dat`). Les potentiels élémentaires λ sont déduits du gaz, et
   aucune phase condensée absente ne doit avoir une force motrice
   g°/RT − Σ a·λ < 0.
3. Si `bprime` bloque ou si le résultat est instable, il énumère les
   assemblages de phases possibles (un polymorphe valide par famille : C,
   Si, SiC, SiO2), avec `bprime` sur un mélange restreint.
4. Là où plusieurs solides coexistent (SiC + C(gr), SiC + SiO2(L) +
   Si(L)…), le solveur échoue même sur le mélange restreint. L'équilibre à
   assemblage fixé est alors résolu directement en Python, en potentiels
   élémentaires, avec les mêmes données, les mêmes masses atomiques et le
   même excès de char que `bprime`. Sur 300 points tirés au hasard, ce
   solveur reproduit `bprime` à 8e-6 près (B'c) et 4e-6 près (h_w).

Bilan sur la grille : 4387 points directs, 289 par énumération, 49 à
assemblage fixé. **Tous sont des équilibres strictement stables.**
`Si3N4(cr)` est exclu du mélange : il bloque systématiquement le solveur, et
la nitruration du SiC sous air est cinétiquement inhibée.

La bibliothèque Mutation++ n'est pas modifiée.

## 3. Régimes (1 atm)

| Tw [K] | phases condensées | mécanisme | B'c |
|--|--|--|--|
| 300 – 1625 | SiC + SiO2 + C(gr) | oxydation passive : la paroi fixe l'oxygène | 0 (prise de masse écrêtée) |
| 1650 – 1950 | SiC + SiO2 | transition passive → active | 0 → 0.17 |
| 1975 – 2200 | SiC | oxydation active SiC + O2 → SiO + CO | **0.292** = y_O·M_SiC/M_O2 |
| 2225 – 3000 | SiC + Si(L) | oxydation, le silicium reste liquide | 0.29 → 1.7 |
| 3025 – 3100 | SiC | sublimation (Si, SiC2, Si2C) | 2.3 → 13.5 |
| ≥ 3125 | gaz seul | sublimation totale | 200 (plafond du solveur) |

La pression décale toutes les transitions. À 0.001 atm, l'oxydation est
active dès 1500 K et le SiC se décompose en laissant un résidu de carbone
dès 2350 K. À 1000 atm, la silice liquide protège jusqu'à ~2700 K.

## 4. Domaine de validité

La table suppose l'équilibre thermochimique à la paroi et une ablation
stationnaire. Pour le SiC, cela ne vaut pas partout :

- **Oxydation active, sublimation** : réactions rapides, limitées par la
  diffusion dans la couche limite. La table B' s'applique comme pour le carbone.
- **Oxydation passive** : B'c = 0 est juste qualitativement (récession
  négligeable). En réalité, la croissance de la silice est limitée par la
  diffusion de l'oxygène à travers la couche d'oxyde. La prise de masse et
  l'épaisseur de SiO2 relèvent d'un modèle cinétique.
- **Transition passive → active** : elle est gouvernée par la cinétique et
  le transport de SiO et de O2 (type Wagner). L'équilibre ne la place pas
  au bon endroit, et le seuil doit venir d'un critère cinétique ou d'essais.
- **Solide résiduel** (Si liquide, carbone) : B'c reste la masse gazéifiée,
  mais la récession du SiC n'en découle plus directement.

La colonne `phases_condensees` de `sic_bprime_phases.csv` indique le régime
de chaque point.

Mise en données détaillée : onglet `SiC` de
[`../mise_en_donnees_xlsx/mise_en_donnees_materiaux.xlsx`](../mise_en_donnees_xlsx/).
