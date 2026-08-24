# Mises en données de pyrolyse du CPh70 : formes en `ρ_i,c`, en `F_i`, en `Δρ_i`

Pendant de `../tacot_bprime/mise_en_donnees_pyrolyse.md`. Le CPh70 emploie
**les mêmes fibres et la même résine que le TACOT** — seules les proportions et
la porosité changent :

| | TACOT | **CPh70** |
|---|---|---|
| fibres / résine (volume solide) | 50 / 50 | **70 / 30** |
| porosité | 0.80 | **0.01** |
| ε_f | 0.1000 | **0.6930** |
| ε_m | 0.1000 | **0.2970** |

```bash
python tacot_bprime/tacot_mise_en_donnees.py CPh70
```

---

## 0. Ce qui change, ce qui ne change pas

**La cinétique ne change pas du tout.** `A`, `E/R`, `ψ`, `T_reac` sont ceux de
Goldstein pour la résine du TACOT (`TACOT_3.0.xls`, `Pyrolysis model!A11:G20`).
Seul **`ε_m` change**, et il ne fait que mettre à l'échelle les densités.

**L'avertissement du TACOT s'applique intégralement.** La phase B a une densité
résiduelle `ρ_c` **non nulle** (900 → 600 kg/m³ de matrice). Les formes « en
`F` » et « en `Δρ` » supposent implicitement `ρ_i,c = 0` et exigent donc de
diviser `A_B` par **9**. L'oubli coûte ici **35.4 kg/m³**, soit **20 % du gaz
total libéré** — même proportion que pour le TACOT, puisque c'est la même
résine.

---

## 1. JEU 0 — la forme native

```
d(ρ_i)/dt = − A_i · ρ_i,v · [(ρ_i − ρ_i,c)/ρ_i,v]^ψ_i · exp(−E_i/(R·T))
            pour T > T_reac, 0 sinon
```

### (a) Base résine — identique au TACOT

| phase | ρ_i,v | ρ_i,c | A [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|---|
| A | 300.0 | 0.0 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 900.0 | **600.0** | 4.4800e+09 | 20444.444444 | 3.0 | 555.5556 |
| C | 1600.0 | 1600.0 | 0 | — | — | — |

`Σ ρ_i,v = 1200 kg/m³` = densité intrinsèque de la résine.

### (b) Base composite — multiplier par `ε_m = 0.297`

| phase | ρ_i,v | ρ_i,c | A [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|---|
| A | 89.1 | 0.0 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 267.3 | 178.2 | 4.4800e+09 | 20444.444444 | 3.0 | 555.5556 |
| C (fibres) | **1108.8** | 1108.8 | 0 | — | — | — |

Somme vierge **1465.2**, somme char **1287.0** kg/m³.

> **Les `A` sont inchangés** entre les deux bases : l'équation native est
> invariante par mise à l'échelle commune de `(ρ_i, ρ_i,v, ρ_i,c)`.

---

## 2. La correction qu'exigent les jeux A et B

| phase | (ρ_v−ρ_c)/ρ_v | facteur `[·]^(ψ−1)` | A publié | **A à employer** |
|---|---|---|---|---|
| A | 1.000000 | 1.000000 | 1.2000e+04 | **1.2000e+04** |
| B | 0.333333 | **0.111111** | 4.4800e+09 | **4.9778e+08** |

Identique au TACOT — c'est une propriété de la résine, pas du composite.

---

## 3. JEU A — Arrhenius + `F_i`

```
dx_i/dt = − A_i · x_i^ψ_i · exp(−E_i/(R·T))          x_i(0) = 1
ρ(T)    = ρ_vierge · [ 1 − Σ_i F_i · (1 − x_i) ]
```

| | valeur |
|---|---|
| ρ_vierge | 1465.20 kg/m³ |
| ρ_char | 1287.00 kg/m³ |
| porosité vierge | 0.0100 |

| i | F_i [-] | A_i [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|
| A | 0.060810811 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 0.060810811 | **4.9778e+08** | 20444.444444 | 3.0 | 555.5556 |
| **Σ** | **0.121621622** | | | | |

Contrôle : `1 − ΣF_i = 0.878378378 = ρ_char/ρ_vierge`.

> À comparer au TACOT, dont `ΣF = 0.214286`. Le CPh70 perd **près de deux fois
> moins** de sa masse en pyrolyse : il est riche en fibres, qui ne pyrolysent
> pas.

---

## 4. JEU B — Arrhenius + `Δρ_i`

```
d(ρ_i)/dt = − A_i · Δρ_i · (ρ_i/Δρ_i)^ψ_i · exp(−E_i/(R·T))     ρ_i(0) = Δρ_i
ρ(T) = ρ_inerte + Σ_i ρ_i                        ρ_inerte = ρ_char
```

| | valeur |
|---|---|
| ρ_vierge | 1465.20 kg/m³ |
| ρ_inerte (= ρ_char) | 1287.00 kg/m³ — dont **1108.8** de fibres et **178.2** de char de résine |

| i | Δρ_i [kg/m³] | A_i [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|
| A | 89.100000 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 89.100000 | **4.9778e+08** | 20444.444444 | 3.0 | 555.5556 |
| **Σ** | **178.200000** | | | | |

Contrôle : `1287.00 + 178.20 = 1465.20 = ρ_vierge`.

> **Si le code réclame un couple `(ρ_i,v ; ρ_i,c)`** : soit le **JEU 0** tel
> quel avec les `A` publiés, soit `ρ_i,v = Δρ_i` et `ρ_i,c = 0` **avec les `A`
> corrigés**. Panacher les deux est l'erreur classique, et la phase A étant
> juste dans les deux lectures, elle passe inaperçue.

---

## 5. Vérification

Trois intégrations **indépendantes**, ATG à 20 K/min, densité du composite en
kg/m³ :

| T [°C] | JEU 0 natif | JEU A | JEU B | *sans la correction* |
|---|---|---|---|---|
| 300 | 1442.3792 | 1442.3792 | 1442.3792 | *1442.3756* |
| 400 | 1402.7394 | 1402.7394 | 1402.7394 | *1401.3516* |
| 500 | 1376.9395 | 1376.9395 | 1376.9395 | ***1346.0847*** |
| 550 | 1350.6300 | 1350.6300 | 1350.6300 | ***1317.6825*** |
| 600 | 1323.6866 | 1323.6866 | 1323.6866 | *1303.2941* |
| 700 | 1298.8925 | 1298.8925 | 1298.8925 | *1292.8538* |
| 900 | 1289.3948 | 1289.3948 | 1289.3948 | *1288.5460* |

Écart maximal au jeu natif : **4.8e-12** (jeu A) et **5.0e-12** kg/m³ (jeu B).

Sans la correction : **35.4 kg/m³**, soit **20 %** du gaz total libéré.

---

## 6. Ce qui accompagne le bloc cinétique

```
gaz de pyrolyse    C:0.206, H:0.679, O:0.115   (fractions molaires élémentaires)
char               C:1.0
bord de couche l.  N:0.79, O:0.21              (air)
ρ_v / ρ_c          1465.20 / 1287.00 kg/m³
ε_f / ε_m          0.6930 / 0.2970     porosité 0.0100
couplage k         0.1385      = B'g/B'c en ablation stationnaire
```

**Le gaz de pyrolyse et le char sont ceux du TACOT** — c'est la même résine et
les mêmes fibres. C'est pourquoi `data/mixtures/cph70-air_25.xml` est identique
à `tacot-air_25.xml`, et pourquoi **la table B' est rigoureusement la même**
(écart 0.000e+00, cf. `../tacot_bprime/autre_materiau.md`).

---

## 7. Ce que le CPh70 change vraiment

| | TACOT | CPh70 | rapport |
|---|---|---|---|
| ρ vierge | 280.0 | **1465.2** | ×5.2 |
| ρ char | 220.0 | **1287.0** | ×5.9 |
| Σ F_i (perte de masse) | 0.214286 | **0.121622** | ×0.57 |
| gaz libéré | 60.0 kg/m³ | **178.2 kg/m³** | ×2.97 |
| couplage k = B'g/B'c | 0.2727 | **0.1385** | ×0.51 |
| table B' | — | **identique** | ×1 |

Trois lectures de ce tableau :

- **Le CPh70 libère trois fois plus de gaz par m³**, mais **deux fois moins en
  proportion de sa masse** — il est dense et riche en fibres.
- **Le couplage `k` est divisé par deux** : à `B'c` donné, le CPh70 souffle
  deux fois moins de gaz de pyrolyse. Son point de fonctionnement sur la table
  B' est donc différent, même si la table ne l'est pas.
- **À `B'c` égal, il recule 5.85 fois moins vite** (1287/220), la récession
  étant `ṡ = B'c·ṁ_e/ρ_c`.

---

## 8. Fichiers

| Fichier | Rôle |
|---|---|
| `../tacot_bprime/tacot_mise_en_donnees.py` | émet les trois jeux (`TACOT` ou `CPh70`) et les vérifie |
| `../tacot_bprime/mise_en_donnees_pyrolyse.md` | le même document pour le TACOT |
| `../tacot_bprime/autre_materiau.md` | pourquoi la table B' ne change pas |
| `../tacot_bprime/material_response.py` | bilan de phase solide, lignes de fonctionnement |
| `README.md` | table B' et validation du CPh70 |
| `../data/mixtures/cph70-air_25.xml` | mise en données du mélange (identique au TACOT) |
