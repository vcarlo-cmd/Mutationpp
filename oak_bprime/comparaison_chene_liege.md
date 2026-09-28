# Chêne vs liège/phénolique P50 — comparaison des tables B'

Comparaison de la table B' du bois de chêne (`oak_bprime/`) avec celle du
liège/phénolique **P50 telle que déjà calculée** (`../cork_bprime/`).
Hypothèse commune : **rendement en char 20 %** pour les deux matériaux.

Tous les chiffres ci-dessous sont relus dans les CSV versionnés par
`comparaison_chene_liege.py`, sans nouveau calcul d'équilibre :

```bash
python comparaison_chene_liege.py   # tableaux + comparaison_chene_liege.{csv,png}
```

---

## 1. Même rendement en char, même k

| | chêne | liège/phénolique P50 |
|---|---|---|
| rendement en char | **20 %** (donnée de l'énoncé) | **20 %** (TGA P50, argon) |
| bilan sur 100 g vierge | 20 g char + 80 g gaz | 20 g char + 80 g gaz |
| k = B'g/B'c = (1−y)/y | **4,0** | **4,0** |
| char | C:1.0 | C:1.0 |
| bord de couche limite | air N:0,79 O:0,21 | air N:0,79 O:0,21 |
| mélange paroi | 25 espèces + C(gr) | même liste d'espèces |
| ρ vierge | 700 kg/m³ (ordre de grandeur) | 465,6 kg/m³ (mesurée) |

- **Char et air identiques** : à B'g = 0 les deux tables **coïncident
  exactement** (B'c = 0,0874 à 300 K, limite C + O₂ → CO₂ ; plateau à 0,175
  de 1500 à 2500 K ; sublimation vers 3500-4000 K).
- **Même k = 4** : pour chaque gramme de char consommé en paroi, les deux
  matériaux soufflent 4 g de gaz de pyrolyse.

À rendement en char égal, **la seule différence entre les deux tables est la
composition du gaz de pyrolyse**.

## 2. La chimie du gaz de pyrolyse

| gaz de pyrolyse | C | H | O | O/C |
|---|---|---|---|---|
| chêne (`oak_pyro`) | 0,221 | 0,536 | 0,243 | **1,10** |
| liège/phénolique P50 (`cork_pyro`) | 0,287 | 0,592 | 0,121 | **0,42** |

- **Gaz du liège P50 pauvre en oxygène** (O/C = 0,42) : il ne peut pas
  oxyder le char et dilue l'oxygène de l'air. Le soufflage **protège** le
  char.
- **Gaz du chêne riche en oxygène** (O/C = 1,10, H₂O et CO₂ dominants) : au-
  dessus de ~1100 K il **oxyde** le char. Le soufflage **augmente** B'c.

B'c à 1 atm, B'g imposé (chêne / liège P50) :

| T [K] | B'g = 0 | B'g = 0,5 | B'g = 2 | B'g = 5 |
|---|---|---|---|---|
| 1000 | 0,154 / 0,154 | 0,128 / 0,000 | 0,056 / 0 | 0 / 0 |
| 1500 | 0,175 / 0,175 | 0,194 / 0,009 | 0,250 / 0 | 0,363 / 0 |
| 2000 | 0,175 / 0,175 | 0,199 / 0,014 | 0,260 / 0 | 0,380 / 0 |
| 2500 | 0,175 / 0,175 | 0,220 / 0,040 | 0,313 / 0 | 0,478 / 0 |
| 3000 | 0,177 / 0,177 | 0,274 / 0,105 | 0,463 / 0 | 0,791 / 0 |

Pour le liège P50, B'c est nul dès B'g = 2 jusqu'à ~3300 K ; pour le chêne,
il est multiplié par 2 à 4,5 à B'g = 5 selon la température.

## 3. Point de fonctionnement stationnaire (B'g = 4·B'c, 1 atm)

| T [K] | B'c chêne | B'c liège P50 | B'g chêne | B'g liège P50 | rapport B'c |
|---|---|---|---|---|---|
| 1000 | 0,127 | 0,060 | 0,51 | 0,24 | ×2,1 |
| 1500 | 0,206 | 0,075 | 0,82 | 0,30 | ×2,7 |
| 2000 | 0,214 | 0,077 | 0,85 | 0,31 | ×2,8 |
| 2500 | 0,254 | 0,086 | 1,02 | 0,34 | ×3,0 |
| 3000 | 0,431 | 0,114 | 1,72 | 0,46 | ×3,8 |

Référence sans pyrolyse (B'g = 0) : 0,154 / 0,175 / 0,177 à 1000 / 2000 /
3000 K. Le soufflage du liège P50 **divise B'c par 2,6 à 1000 K, 2,3 à
2000 K et 1,6 à 3000 K** ; celui du chêne le réduit encore un peu à 1000 K
(−17 %) mais **l'augmente de 22 % à 2000 K et le multiplie par 2,4 à
3000 K**.

**Perte de masse.** k étant le même, la perte de masse totale
ṁ/ṁe = (1+k)·B'c est dans le même rapport que B'c : le chêne perd
**2,8 fois plus de masse à 2000 K et 3,8 fois plus à 3000 K**. Ce rapport ne
dépend d'aucune masse volumique.

**Récession.** En régime stationnaire, le front recule de
ṡ/ṁe = (1+k)·B'c/ρ_v :

| T [K] | chêne (ρ_v = 700) | liège P50 (ρ_v = 465,6) | rapport |
|---|---|---|---|
| 1000 | 0,91 | 0,64 | ×1,4 |
| 2000 | 1,53 | 0,83 | ×1,8 |
| 3000 | 3,08 | 1,22 | ×2,5 |

(unités : 10⁻³ m³/kg). La masse volumique plus élevée du chêne compense en
partie : il récède **1,8 fois plus vite à 2000 K et 2,5 fois plus vite à
3000 K**. Ce rapport hérite de l'incertitude sur ρ_v du chêne (700 kg/m³,
ordre de grandeur) ; celui sur la perte de masse, non.

> La colonne `recession_over_mdote_m3_per_kg` des CSV de chaque matériau
> vaut B'c/ρ_c (vitesse de consommation du seul char, avec ρ_c = 255 kg/m³
> supposée pour le chêne et 289,1 kg/m³ mesurée pour le P50). Ce n'est pas
> la grandeur utilisée ici.

Au-delà de ~3200 K à 1 atm, le point stationnaire du chêne sort de la table
(B'g_ss atteint 10 à 3400 K) : la comparaison s'arrête à 3000 K.

## 4. Conclusion

À rendement en char égal (20 %) et donc à soufflage égal (k = 4), le liège/
phénolique P50 ablate **2 à 4 fois moins** (en masse) que le chêne, entre
1000 et 3000 K à 1 atm. Ce n'est pas une question de quantité de char ou de
gaz, mais de **chimie du gaz** : celui du P50 (O/C = 0,42) protège le char,
celui du chêne (O/C = 1,10) l'oxyde.

---

## Sources

- `comparaison_chene_liege.py` → `comparaison_chene_liege.csv` (point
  stationnaire des deux matériaux) et `comparaison_chene_liege.png`
  (B'c à B'g imposé, et point stationnaire, 1 atm)
- `oak_bprime_bc_table.csv`, `oak_bprime_steady_state.csv` (`oak_bprime.py`)
- `../cork_bprime/cork_bprime_bc_table.csv`,
  `../cork_bprime/cork_bprime_steady_state.csv` (`cork_bprime.py`, P50)
