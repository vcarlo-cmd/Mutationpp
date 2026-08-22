# Table B' et enthalpie du gaz de pyrolyse — Zuram sous air

**Zuram** — composite carbone/phénolique VKI, cas de référence AblaNTIS
`carbonPhenolInAir`.

| | |
|--|--|
| mélange paroi | `data/mixtures/zuram-air.xml` (20 espèces gazeuses + C(gr)) |
| mélange gaz de pyrolyse | `data/mixtures/zuram-pyrogas.xml` (mêmes espèces, sans C(gr)) |
| gaz de pyrolyse (`VKIZuramPyroGas`) | C:0.171, H:0.722, O:0.107 ([traçabilité](resine_zuram.md)) |
| arête de couche limite (`BLedge`) | O:0.21, N:0.79 (air) |
| char (`Char`) | C:1.0 |
| table de référence | `Bprime_carbonPhenolInAir_AblaNTIS.txt` (15600 points) |

```bash
cd zuram_bprime
python zuram_bprime.py          # table B'c / h_w + figures
python zuram_validation.py      # comparaison à la table AblaNTIS
python zuram_pyrolysis_gas.py   # h_g, M, Cp, gamma, rho, mu du gaz de pyrolyse
python zuram_variantes.py       # variantes XX/YY : densités, ATG, porosité
python zuram_variantes_bprime.py  # ... et leur vérification par bprime
python zuram_cinetique.py       # transposition des paramètres d'Arrhenius
python zuram_mise_en_donnees.py # les deux jeux de cinétique, base composite
```

> Traçabilité de la résine et de la nomenclature « 18/50 » :
> [`resine_zuram.md`](resine_zuram.md).
> Variantes hypothétiques (14/40, 18/80…) :
> [`variantes_zuram.md`](variantes_zuram.md).

Les trois scripts cherchent `bprime` / `mppequil` dans le `PATH` puis dans
`build/src/apps/`, et positionnent `MPP_DATA_DIRECTORY` sur `data/` s'il n'est
pas déjà défini.

---

## 1. `zuram_bprime.py` — table B' et h_w

Balaye 25 isobares log-espacées de 10⁻³ à 10³ atm, 189 températures de 300 à
5000 K (pas 25 K), pour six valeurs de B'g : 0, 0.1, 0.2, 0.5, 1.0, 2.0.

Chaque CSV contient, en plus de B'c et h_w, les 21 fractions molaires de paroi
et une colonne `regime`.

> **Plafond numérique.** `Thermodynamics::surfaceMassBalance` ajoute une
> quantité **finie** de char, `max(100·B'g, 200)`. Au-delà de la sublimation
> complète, tout ce char passe en phase gazeuse et B'c sature exactement sur
> cette valeur. Ce n'est pas une solution physique mais la signature de
> B'c → ∞. Ces points portent `regime = SATURE` dans les CSV et sont masqués
> sur les figures — les filtrer avec `regime == "OK"`.

## 2. `zuram_validation.py` — validation contre AblaNTIS

Refait le calcul sur la grille **exacte** de la référence
(10 pressions × 24 B'g × 65 températures = 15600 points) et compare point par
point B'c et h_w.

La référence est un fichier texte à 5 colonnes, sans en-tête :

```
P[Pa]   B'g   B'c   T[K]   h_w[J/kg]
```

Les points sont classés en trois régimes : `ACTIF` (comparable), `ZERO`
(la référence impose B'c = 0) et `SATURE` (plafond de char, accord trivial —
isolé pour ne pas flatter les statistiques).

### Résultat

| grandeur | n | écart moyen | écart max | points identiques |
|---|---|---|---|---|
| B'c (régime ACTIF) | 9745 | 0.00e+00 % | **0.00e+00 %** | 100.00 % |
| h_w (régime ACTIF) | 9745 | 3.85e-16 % | **1.46e-14 %** | 96.26 % |
| h_w (ZERO + ACTIF) | 15264 | 3.88e-16 % | 1.46e-14 % | 96.31 % |

- **B'c : identique à tous les chiffres publiés** sur les 9745 points actifs.
- **B'c = 0** retrouvé sur 5519/5519 points du régime ZERO.
- **Plafond de saturation** retrouvé sur 336/336 points.

L'écart résiduel sur h_w (1e-14 %) est du bruit de conversion décimale, pas un
écart physique. La référence n'étant stockée qu'avec 6 chiffres significatifs,
la résolution de la comparaison est de ~1e-4 % ; « écart nul » signifie donc
« identique à tous les chiffres publiés », pas « identique au bit près ».

Cet accord exact était attendu : la table AblaNTIS a été produite avec le même
solveur d'équilibre multiphase et la même liste d'espèces. La validation
confirme que la **mise en données** de `zuram-air.xml` — liste d'espèces,
compositions élémentaires, présence de C(gr) — reproduit bien le cas de
référence.

## 3. `zuram_pyrolysis_gas.py` — enthalpie du gaz de pyrolyse

`h_g` est la grandeur qui manque à la table B' pour fermer le bilan d'énergie
de surface (avec `h_w`) et le terme source de pyrolyse en profondeur
(`h_g - h_s`).

Le gaz de pyrolyse **pur**, à l'équilibre chimique, **sans phase condensée** :
ni air, ni char, ni C(gr) — l'état du gaz avant qu'il n'atteigne la paroi.
25 isobares de 10⁻³ à 10³ atm, T de 200 à 4000 K (pas 25 K).

Quelques valeurs à 1 atm :

| T [K] | h_g [kJ/kg] | M [kg/kmol] | Cp [kJ/kg/K] |
|---|---|---|---|
| 300 | −7787.2 | 19.203 | 1.840 |
| 1000 | −1869.6 | 11.459 | 10.251 |
| 2000 | 5036.2 | 9.617 | 4.453 |
| 3000 | 12392.0 | 8.997 | 14.317 |
| 4000 | 40688.8 | 6.181 | 34.232 |

---

## Fichiers produits

| Fichier | Contenu |
|---------|---------|
| `zuram_bprime_Bg{0p0,0p1,0p2,0p5,1p0,2p0}.csv` | table B' complète : 25 isobares × 189 T × (B'c, h_w, 21 fractions molaires) |
| `zuram_bprime_Bg*.png` | isobares B'c et h_w pour chaque B'g |
| `zuram_bprime_bg_comparison.png` | influence de B'g à 1 atm |
| `zuram_bprime_vs_ref_zuram-air.csv` | comparaison point par point à la référence (15600 lignes) |
| `zuram_bprime_vs_ref_zuram-air.png` | superposition calcul / référence, B'c et h_w |
| `zuram_bprime_err_zuram-air.png` | carte de l'écart relatif sur B'c |
| `zuram_pyrolysis_gas.csv` | h_g, M, Cp, γ, ρ, μ du gaz de pyrolyse (25 isobares × 153 T) |
| `zuram_pyrolysis_gas.png` | les six propriétés par isobare |
| `zuram_pyrolysis_gas_enthalpy.png` | zoom sur h_g et sa sensibilité à la pression |

Voir `tacot_bprime/` pour le même travail sur le TACOT et `cph70_bprime/` pour
le CPh70.
