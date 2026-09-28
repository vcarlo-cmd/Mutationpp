# Chêne / liège seul / liège-phénolique P50 — comparaison des tables B'

Comparaison de la table B' du bois de chêne (`oak_bprime/`) avec celle du
liège/phénolique **P50 telle que déjà calculée** (`../cork_bprime/`) et celle
du **liège seul**, sans résine (`../cork_bprime/cork_seul_bprime.py`).
Hypothèse commune : **rendement en char 20 %** pour les trois matériaux.

Tous les chiffres ci-dessous sont relus dans les CSV versionnés par
`comparaison_chene_liege.py`, sans nouveau calcul d'équilibre :

```bash
python ../cork_bprime/cork_seul_bprime.py   # table du liège seul (bprime requis)
python comparaison_chene_liege.py           # tableaux + comparaison_chene_liege.{csv,png}
```

---

## 1. Même rendement en char, même k

| | chêne | liège seul | liège/phénolique P50 |
|---|---|---|---|
| rendement en char | **20 %** (donnée de l'énoncé) | **20 %** (imposé) | **20 %** (TGA P50, argon) |
| bilan sur 100 g vierge | 20 g char + 80 g gaz | idem | idem |
| k = B'g/B'c = (1−y)/y | **4,0** | **4,0** | **4,0** |
| char | C:1.0 | C:1.0 | C:1.0 |
| bord de couche limite | air N:0,79 O:0,21 | idem | idem |
| mélange paroi | `oak-air` | `corkpure-air` | `cork-air` |
| ρ vierge [kg/m³] | 700 (ordre de grandeur) | 372 (**hypothèse**, voir §3) | 465,6 (mesurée) |

- **Char et air identiques** : à B'g = 0 les trois tables **coïncident
  exactement** (B'c = 0,0874 à 300 K, limite C + O₂ → CO₂ ; plateau à 0,175
  de 1500 à 2500 K ; sublimation vers 3500-4000 K).
- **Même k = 4** : pour chaque gramme de char consommé en paroi, les trois
  matériaux soufflent 4 g de gaz de pyrolyse.

À rendement en char égal, **la seule différence entre les tables est la
composition du gaz de pyrolyse**.

Le liège seul reprend l'analyse élémentaire du liège du P50 (C 62,4 / H 8,5 /
O 28,4 % masse) ; seul le rendement en char est imposé à 20 %. Dans le P50,
le rendement du liège déduit de la TGA est 12,5 % : le cas « liège seul à
20 % » est un cas de comparaison, pas un matériau mesuré.

## 2. La chimie du gaz de pyrolyse

| gaz de pyrolyse | C | H | O | O/C |
|---|---|---|---|---|
| chêne (`oak_pyro`) | 0,221 | 0,536 | 0,243 | **1,10** |
| liège seul (`corkpure_pyro`) | 0,258 | 0,613 | 0,129 | **0,50** |
| liège/phénolique P50 (`cork_pyro`) | 0,287 | 0,592 | 0,121 | **0,42** |

- **Gaz des deux lièges pauvres en oxygène** (O/C ≤ 0,5) : ils ne peuvent pas
  oxyder le char et diluent l'oxygène de l'air. Le soufflage **protège** le
  char.
- **La résine abaisse O/C de 0,50 à 0,42** : le gaz de la novolac (C₇H₆O,
  O/C = 0,39 après carbonisation) est plus carboné que celui du liège. Le P50
  protège donc un peu mieux que le liège seul.
- **Gaz du chêne riche en oxygène** (O/C = 1,10, H₂O et CO₂ dominants) : au-
  dessus de ~1100 K il **oxyde** le char. Le soufflage **augmente** B'c.

B'c à 1 atm, B'g imposé (chêne / liège seul / liège P50) :

| T [K] | B'g = 0 (les trois) | B'g = 0,5 | B'g = 2 | B'g = 5 |
|---|---|---|---|---|
| 1000 | 0,154 | 0,128 / 0,000 / 0,000 | 0,056 / 0 / 0 | 0 / 0 / 0 |
| 1500 | 0,175 | 0,194 / 0,042 / 0,009 | 0,250 / 0 / 0 | 0,363 / 0 / 0 |
| 2000 | 0,175 | 0,199 / 0,047 / 0,014 | 0,260 / 0 / 0 | 0,380 / 0 / 0 |
| 2500 | 0,175 | 0,220 / 0,074 / 0,040 | 0,313 / 0 / 0 | 0,478 / 0 / 0 |
| 3000 | 0,177 | 0,274 / 0,143 / 0,105 | 0,463 / 0 / 0 | 0,791 / 0 / 0 |

À B'g = 2, B'c est nul jusqu'à ~3100 K pour le liège seul et ~3300 K pour le
P50 ; pour le chêne, il est multiplié par 2 à 4,5 à B'g = 5 selon la
température. L'écart entre les deux lièges n'apparaît qu'à faible soufflage
(B'g = 0,5).

## 3. Point de fonctionnement stationnaire (B'g = 4·B'c, 1 atm)

| T [K] | B'c chêne | B'c liège seul | B'c liège P50 | B'c(B'g = 0) |
|---|---|---|---|---|
| 1000 | 0,127 | 0,066 | 0,060 | 0,154 |
| 1500 | 0,206 | 0,085 | 0,075 | 0,175 |
| 2000 | 0,214 | 0,087 | 0,077 | 0,175 |
| 2500 | 0,254 | 0,098 | 0,086 | 0,175 |
| 3000 | 0,431 | 0,137 | 0,114 | 0,177 |

B'g stationnaire = 4·B'c : 0,85 / 0,35 / 0,31 à 2000 K, 1,72 / 0,55 / 0,46 à
3000 K (chêne / liège seul / P50).

Le soufflage des deux lièges **divise B'c par 1,8 à 2,6** jusqu'à 2500 K (par
1,3 pour le liège seul et 1,6 pour le P50 à 3000 K) ; celui du chêne le réduit un peu à 1000 K (−17 %) mais
**l'augmente de 22 % à 2000 K et le multiplie par 2,4 à 3000 K**.

**Perte de masse.** k étant le même, la perte de masse totale
ṁ/ṁe = (1+k)·B'c est dans le même rapport que B'c. Ce rapport ne dépend
d'aucune masse volumique :

| T [K] | chêne / P50 | liège seul / P50 | chêne / liège seul |
|---|---|---|---|
| 1000 | ×2,1 | ×1,11 | ×1,9 |
| 1500 | ×2,7 | ×1,13 | ×2,4 |
| 2000 | ×2,8 | ×1,13 | ×2,5 |
| 2500 | ×3,0 | ×1,15 | ×2,6 |
| 3000 | ×3,8 | ×1,20 | ×3,1 |

**Récession.** En régime stationnaire, le front recule de
ṡ/ṁe = (1+k)·B'c/ρ_v (10⁻³ m³/kg) :

| T [K] | chêne (ρ_v = 700) | liège seul (ρ_v = 372) | liège P50 (ρ_v = 465,6) |
|---|---|---|---|
| 1000 | 0,91 | 0,89 | 0,64 |
| 2000 | 1,53 | 1,17 | 0,83 |
| 3000 | 3,08 | 1,84 | 1,22 |

- ρ_v du chêne (700 kg/m³) est un ordre de grandeur.
- ρ_v du liège seul n'est pas mesurée : 372 kg/m³ = 0,8 × 465,6, c'est-à-
  dire le P50 dont on retire la résine à volume constant. Ce n'est qu'une
  hypothèse de travail : un aggloméré de liège sans liant peut être nettement
  plus léger.
- Seule la ligne P50 repose sur une masse volumique mesurée ; les rapports de
  perte de masse ci-dessus, eux, sont exacts quelle que soit ρ_v.

> La colonne `recession_over_mdote_m3_per_kg` des CSV du chêne et du P50
> vaut B'c/ρ_c (vitesse de consommation du seul char, avec ρ_c = 255 kg/m³
> supposée pour le chêne et 289,1 kg/m³ mesurée pour le P50). Ce n'est pas
> la grandeur utilisée ici. Le CSV du liège seul donne directement
> (1+k)·B'c (colonne `mass_loss_over_mdote`).

Au-delà de ~3200 K à 1 atm, le point stationnaire du chêne sort de la table
(B'g_ss atteint 10 à 3400 K) : la comparaison s'arrête à 3000 K.

## 4. Conclusion

À rendement en char égal (20 %), donc à soufflage égal (k = 4), entre 1000 et
3000 K à 1 atm :

1. **Chêne vs liège** : le chêne perd **2 à 4 fois plus de masse** que les
   deux lièges. C'est la chimie du gaz qui décide : celui du chêne
   (O/C = 1,10) oxyde le char, ceux du liège (O/C ≤ 0,5) le protègent.
2. **Liège seul vs P50** : à rendement égal, la résine ne gagne que **11 à
   20 %** de perte de masse, par un gaz un peu plus carboné (O/C 0,42 au lieu
   de 0,50). L'effet chimique de la résine est donc secondaire ; son apport
   principal est ailleurs : elle relève le rendement en char réel du
   composite (12,5 % pour le liège dans le P50 → 20 % pour le composite) et
   donne un char cohérent, deux effets que ce cas à rendement imposé ne
   représente pas.

---

## Sources

- `comparaison_chene_liege.py` → `comparaison_chene_liege.csv` (point
  stationnaire des trois matériaux) et `comparaison_chene_liege.png`
  (B'c à B'g imposé, et point stationnaire, 1 atm)
- `oak_bprime_bc_table.csv`, `oak_bprime_steady_state.csv` (`oak_bprime.py`)
- `../cork_bprime/cork_bprime_bc_table.csv`,
  `../cork_bprime/cork_bprime_steady_state.csv` (`cork_bprime.py`, P50)
- `../cork_bprime/cork_seul_bprime_bc_table.csv`,
  `../cork_bprime/cork_seul_bprime_steady_state.csv`
  (`cork_seul_bprime.py`, mélange `data/mixtures/corkpure-air.xml`)
