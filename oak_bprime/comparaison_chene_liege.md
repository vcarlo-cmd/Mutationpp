# Chêne vs liège/phénolique — comparaison des tables B'

Comparaison de la table B' du bois de chêne (`oak_bprime/`) à celle du
liège/phénolique 80/20 (`../cork_bprime/`), et question annexe : quel est
l'intérêt d'ajouter ~20 % de résine phénolique au liège ?

---

## 1. Ce que les deux matériaux ont en commun

| | chêne | liège/phénolique P50 |
|---|---|---|
| rendement en char | 20 % | 20 % |
| k = B'g/B'c = (1−y)/y | **4,0** | **4,0** |
| char | C:1.0 | C:1.0 |
| mélange paroi | 25 espèces + C(gr) | même liste d'espèces |

- **Char et air identiques** ⇒ à B'g = 0 les deux tables B' **coïncident
  exactement** : B'c = 0,0874 à 300 K (limite C + O₂ → CO₂), plateau à 0,175
  de 1500 à 2500 K, sublimation vers 3500-4000 K.
- **Même k = 4** : les deux matériaux gardent 20 % de leur masse en char et
  produisent 4 g de gaz par gramme de char.

La seule différence entre les deux tables vient donc de la **composition du
gaz de pyrolyse**, pas de la quantité de gaz produite.

## 2. Ce qui diffère : la chimie du gaz de pyrolyse

| Gaz de pyrolyse | C | H | O | O/C |
|---|---|---|---|---|
| liège/phénolique P50 (`cork_pyro`) | 0,287 | 0,592 | 0,121 | **0,42** |
| chêne (`oak_pyro`) | 0,221 | 0,536 | 0,243 | **1,10** |

- **Le gaz du liège manque d'oxygène** : il ne peut pas oxyder le char et
  dilue seulement l'oxygène de l'air. Le soufflage **protège** le char : à
  1 atm, B'c tombe à 0 dès B'g ≈ 2, et ce jusqu'à 3000 K.
- **Le gaz du chêne a de l'oxygène en excès** (H₂O, CO₂ dominants) : il
  **oxyde** le char au-dessus de ~1100 K. Le soufflage **augmente** B'c au
  lieu de le faire baisser — comportement inverse du liège.

Effet à 1 atm (voir `oak_bprime_bg_comparison.png` / `cork_bprime_bg0p5.png`,
et la figure combinée `chene_vs_liege_bprime.png`) :

| T [K] | B'c(0) | B'c(B'g=5), chêne | B'c(B'g=5), liège |
|---|---|---|---|
| 1000 | 0,154 | 0 | 0 |
| 2000 | 0,175 | 0,380 | 0 |
| 3000 | 0,177 | 0,791 | 0,020 |

## 3. Point de fonctionnement stationnaire (1 atm)

| T [K] | B'c chêne | B'c liège P50 | rapport |
|---|---|---|---|
| 1000 | 0,127 | 0,060 | ×2,1 |
| 2000 | 0,213 | 0,077 | ×2,8 |
| 3000 | 0,430 | 0,114 | ×3,8 |

Récession estimée ṡ/ṁe = (1+k)·B'c/ρ_v (ρ_v = 700 kg/m³ pour le chêne,
465,6 kg/m³ pour le liège P50, mesurée) : le chêne récède environ **1,8 fois
plus vite à 2000 K et 2,5 fois plus vite à 3000 K**, malgré une masse
volumique plus élevée. Pour la même quantité de gaz produite, le liège
ablate nettement mieux, parce que la chimie de son gaz protège le char au
lieu de l'attaquer.

---

## 4. Intérêt de la résine (~20 %) dans le liège/phénolique

Comparaison, à rendement en char du liège fixé (12,5 %), entre le liège seul
et le composite 80 % liège / 20 % résine (P50) :

| | liège seul | liège + 20 % résine (P50) |
|---|---|---|
| gaz de pyrolyse | C:0,290 H:0,587 O:0,124 | C:0,287 H:0,592 O:0,121 |
| rendement en char du composite | 12,5 % | **20 %** |
| k = B'g/B'c | **7** | 4 |
| B'c stationnaire, 1 atm, 2000 / 3000 K | 0,054 / 0,086 | 0,077 / 0,114 |
| (1+k)·B'c/ρ_v, 2000 / 3000 K [10⁻³ m³/kg] | 1,16 / 1,85 * | **0,83 / 1,22** |

\* ρ_v du liège seul supposée = 0,8 × 465,6 = 372 kg/m³ (hypothèse : le
volume ne change pas quand on retire la résine — pas une mesure).

**La résine ne change presque pas la chimie du gaz** — les deux gaz ont le
même O/C (0,42-0,43) et les tables B'c(T, B'g) coïncident à 0,001 près. Son
effet porte sur la **quantité de char produite** :

1. **Plus de char, moins de récession.** La résine laisse 50 % de char
   contre 12,5 % pour le liège seul. Le rendement composite passe de 12,5 %
   à 20 %, k de 7 à 4 : environ −30 % de récession à 2000 K et −35 % à
   3000 K, malgré un B'c stationnaire plus élevé (moins de soufflage
   protecteur).
2. **Tenue mécanique du char** (hors modèle B'). Le char de liège seul est
   pulvérulent (cellules effondrées). La résine forme un squelette carboné
   continu qui lie les granulés et résiste au cisaillement de l'écoulement —
   sans elle, risque d'érosion mécanique / écaillage non capturé par la
   thermochimie de surface.
3. **Coût.** Matériau plus dense, plus conducteur, soufflage moins fort
   (k = 4 au lieu de 7) : la résine échange un peu d'isolation contre la
   tenue du char.

**Conclusion** : 20 % de résine ne modifie pas la chimie de paroi (même O/C
du gaz), mais transforme le liège — bon isolant qui se désagrège — en
ablateur au char cohérent, avec moins de récession.

---

## Sources

- `oak_bprime/oak_bprime.py`, `oak_bprime/oak_bprime_bc_table.csv`,
  `oak_bprime/oak_bprime_steady_state.csv`
- `../cork_bprime/cork_bprime.py`, `../cork_bprime/cork_bprime_bc_table.csv`,
  `../cork_bprime/cork_bprime_steady_state.csv`
- cas « liège seul » : calcul ad hoc avec `cork_pyrolysis_data.py`
  (`composite_balance(w_cork=1.0, w_resin=0.0)`) et une table B' générée à la
  volée avec le mélange `corkpure-air` (composition C:0,290 H:0,587 O:0,124,
  char C:1.0) — **non versionné dans le dépôt**, reproductible depuis
  `cork_pyrolysis_data.py`.
- figure `chene_vs_liege_bprime.png` — générée pour cette comparaison,
  disponible dans le scratchpad de session (à republier si besoin).
