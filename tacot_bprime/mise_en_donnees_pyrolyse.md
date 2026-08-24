# Mises en données de pyrolyse du TACOT : formes en `ρ_i,c`, en `F_i`, en `Δρ_i`

Pendant de `../zuram_bprime/variantes_zuram.md` § 5. Trois écritures du même
bloc de pyrolyse, toutes vérifiées équivalentes :

| Jeu | Ce qu'il porte | A à employer |
|---|---|---|
| **0 — natif** | `(ρ_i,v ; ρ_i,c)` — ce que publie le classeur | **ceux du classeur** |
| **A** | Arrhenius + `F_i` (fractions, base composite) | **corrigés** |
| **B** | Arrhenius + `Δρ_i` (kg/m³ de composite) | **corrigés** |

Source : `TACOT_3.0.xls`, feuille `Pyrolysis model`, cellules `A11:G20` et
`A31:A34` (modèle de Goldstein, 1965/1969).

```bash
python tacot_bprime/tacot_mise_en_donnees.py
```

---

## 0. L'avertissement, d'abord

**Le TACOT n'est pas le ZURAM.** Sa phase B a une densité résiduelle `ρ_c`
**non nulle** (900 → 600 kg/m³). Les formes « en `F` » et « en `Δρ` », qui
supposent implicitement `ρ_i,c = 0`, ne peuvent donc **pas** reprendre les `A`
du classeur tels quels : il faut diviser `A_B` par **9**.

Oublier cette correction coûte **11.9 kg/m³**, soit **20 % du gaz total
libéré** — une erreur du premier ordre sur le soufflage pyrolytique, pas un
raffinement.

---

## 1. JEU 0 — la forme native, telle que publiée

```
d(ρ_i)/dt = − A_i · ρ_i,v · [(ρ_i − ρ_i,c)/ρ_i,v]^ψ_i · exp(−E_i/(R·T))
            pour T > T_reac, 0 sinon
```

### (a) Base résine — kg/m³ de matrice

| phase | ρ_i,v | ρ_i,c | A [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|---|
| A | 300.0 | 0.0 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 900.0 | **600.0** | 4.4800e+09 | 20444.444444 | 3.0 | 555.5556 |
| C | 1600.0 | 1600.0 | 0 | — | — | — |

`Σ ρ_i,v = 1200 kg/m³` = densité intrinsèque de la résine. La phase C est celle
des fibres : le classeur note explicitement *« The carbon fibers do not
pyrolyse »*.

### (b) Base composite — multiplier par `ε_m = 0.1`

| phase | ρ_i,v | ρ_i,c | A [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|---|
| A | 30.0 | 0.0 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 90.0 | 60.0 | 4.4800e+09 | 20444.444444 | 3.0 | 555.5556 |
| C (fibres) | 160.0 | 160.0 | 0 | — | — | — |

Somme vierge **280.0**, somme char **220.0** — les valeurs de `A32:A33`.

> **Les `A` sont inchangés.** L'équation native est invariante par mise à
> l'échelle commune de `(ρ_i, ρ_i,v, ρ_i,c)` : en posant `ρ' = ε ρ`, les deux
> membres se multiplient par `ε`. Passer de la base résine à la base composite
> ne touche donc pas la cinétique.

---

## 2. La correction qu'exigent les jeux A et B

En posant `y = (ρ_i − ρ_i,c)/(ρ_i,v − ρ_i,c) ∈ [0,1]` dans l'équation native :

```
dy/dt = − A_i · [(ρ_i,v − ρ_i,c)/ρ_i,v]^(ψ−1) · y^ψ · exp(−E_i/(R·T))
```

Le crochet ne vaut 1 que si `ρ_i,c = 0`. Sinon il faut le porter dans `A` :

| phase | (ρ_v−ρ_c)/ρ_v | facteur | A publié | **A à employer** |
|---|---|---|---|---|
| A | 1.000000 | 1.000000 | 1.2000e+04 | **1.2000e+04** |
| B | 0.333333 | **0.111111** | 4.4800e+09 | **4.9778e+08** |

La phase A est intacte — son `ρ_c` était déjà nul. Seule la phase B change,
d'un facteur 9.

---

## 3. JEU A — Arrhenius + `F_i`

```
dx_i/dt = − A_i · x_i^ψ_i · exp(−E_i/(R·T))          x_i(0) = 1
ρ(T)    = ρ_vierge · [ 1 − Σ_i F_i · (1 − x_i) ]
```

| | valeur |
|---|---|
| ρ_vierge | 280.00 kg/m³ |
| ρ_char | 220.00 kg/m³ |
| porosité vierge | 0.8000 |

| i | F_i [-] | A_i [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|
| A | 0.107142857 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 0.107142857 | **4.9778e+08** | 20444.444444 | 3.0 | 555.5556 |
| **Σ** | **0.214285714** | | | | |

Contrôle : `1 − ΣF_i = 0.785714286 = ρ_char/ρ_vierge`.

---

## 4. JEU B — Arrhenius + `Δρ_i`

```
d(ρ_i)/dt = − A_i · Δρ_i · (ρ_i/Δρ_i)^ψ_i · exp(−E_i/(R·T))     ρ_i(0) = Δρ_i
ρ(T) = ρ_inerte + Σ_i ρ_i                        ρ_inerte = ρ_char
```

| | valeur |
|---|---|
| ρ_vierge | 280.00 kg/m³ |
| ρ_inerte (= ρ_char) | 220.00 kg/m³ — dont **160** de fibres et **60** de char de résine |

| i | Δρ_i [kg/m³] | A_i [s⁻¹] | E/R [K] | ψ | T_reac [K] |
|---|---|---|---|---|---|
| A | 30.000000 | 1.2000e+04 | 8555.555556 | 3.0 | 333.3333 |
| B | 30.000000 | **4.9778e+08** | 20444.444444 | 3.0 | 555.5556 |
| **Σ** | **60.000000** | | | | |

Contrôle : `220.00 + 60.00 = 280.00 = ρ_vierge`.

> **Si le code réclame un couple `(ρ_i,v ; ρ_i,c)`** : soit vous employez le
> **JEU 0** tel quel avec les `A` publiés, soit vous posez `ρ_i,v = Δρ_i` et
> `ρ_i,c = 0` **avec les `A` corrigés**. Panacher les deux est l'erreur
> classique — et elle passe inaperçue parce que la phase A, elle, est correcte
> dans les deux lectures.

---

## 5. Vérification

Trois intégrations **indépendantes**, ATG à 20 K/min, densité du composite en
kg/m³ :

| T [°C] | JEU 0 natif | JEU A | JEU B | *sans la correction* |
|---|---|---|---|---|
| 300 | 272.3162 | 272.3162 | 272.3162 | *272.3150* |
| 400 | 258.9695 | 258.9695 | 258.9695 | *258.5022* |
| 500 | 250.2826 | 250.2826 | 250.2826 | ***239.8938*** |
| 550 | 241.4243 | 241.4243 | 241.4243 | ***230.3308*** |
| 600 | 232.3524 | 232.3524 | 232.3524 | *225.4862* |
| 700 | 224.0042 | 224.0042 | 224.0042 | *221.9710* |
| 900 | 220.8063 | 220.8063 | 220.8063 | *220.5205* |

Écart maximal au jeu natif : **1.2e-12** (jeu A) et **8.2e-13** kg/m³ (jeu B).

Sans la correction : **11.9 kg/m³**, soit **20 %** du gaz total libéré, avec un
pic de DTG décalé de 63 °C (549.9 → 486.5 °C).

---

## 6. Ce qui accompagne le bloc cinétique

```
gaz de pyrolyse    C:0.206, H:0.679, O:0.115   (fractions molaires élémentaires)
char               C:1.0
bord de couche l.  N:0.79, O:0.21              (air)
ρ_v / ρ_c          280 / 220 kg/m³
ε_f / ε_m          0.1 / 0.1     porosité 0.80
```

Deux points à ne pas manquer :

- **Le classeur d'origine ne porte aucune enthalpie de pyrolyse.** Les
  −4×10⁶ J/kg qu'on rencontre parfois viennent de la base VKI
  (`Tacot_Zuram_Calcarb_database`, onglet `TACOT_3_0`), pas de `TACOT_3.0.xls`.
- **Deux conventions de pyrolyse coexistent et ne disent pas la même chose sur
  la porosité.** Goldstein décrit une perte de **densité** à volume constant
  (matrice 1200 → 600, `ε_m` reste 0.1, porosité reste 0.80), tandis que
  `A32:A34` décrit une perte de **volume** à densité constante (`ε_m` 0.10 →
  0.05, porosité 0.80 → 0.85). Les deux donnent `ρ_char = 220`, mais pas la
  même perméabilité. Choisir, et s'y tenir.

---

## 7. Comparaison au ZURAM

| | TACOT | ZURAM 18/50 |
|---|---|---|
| réactions | 2 | 4 |
| `ρ_i,c` non nul ? | **oui** (phase B) | non |
| correction sur A | **oui**, facteur 9 sur B | aucune |
| seuils `T_reac` | 333.3 et 555.6 K | aucun |
| Σ F_i (composite) | 0.214286 | 0.190117 |
| rendement en char de la résine | 0.50 | 0.6198 |
| ρ_v / ρ_c | 280 / 220 | 419 / 337.6 (mesuré) |
| gaz libéré | 60 kg/m³ (21.4 %) | 81.4 kg/m³ (19.4 %) |

C'est précisément parce que le ZURAM n'a pas de `ρ_i,c` que sa mise en données
se transpose sans correction — voir
`../zuram_bprime/variantes_zuram.md` § 4, qui utilise justement le TACOT comme
contre-épreuve de la forme en `f`.

---

## 8. Fichiers

| Fichier | Rôle |
|---|---|
| `tacot_mise_en_donnees.py` | émet les trois jeux et vérifie leur équivalence |
| `resine_tacot.md` | traçabilité de la résine jusqu'à Sykes (1967) |
| `mise_en_donnees_xml.md` | ce qu'il faut renseigner dans un XML de mélange |
| `../zuram_bprime/variantes_zuram.md` | le même exercice pour le ZURAM |
