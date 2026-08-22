# Variantes ZURAM XX/YY : 14/40, 18/80 — ce qui change et ce qui ne change pas

Question de conception : en gardant les **mêmes fibres Calcarb** et la **même
résine phénolique**, et en lisant le nom comme au § 1 de `resine_zuram.md` —

```
XX = masse volumique apparente de la préforme, en 0.01 g/cm³
YY = teneur en résine du composite vierge, en % MASSE
```

— que deviennent la masse, la table B', les ATG et la porosité ?

```bash
python zuram_bprime/zuram_variantes.py            # bilan de phase solide
python zuram_bprime/zuram_variantes_bprime.py     # vérification par bprime
```

---

## 0. Réponse courte

| Question | Réponse |
|---|---|
| Plus ou moins dense ? | **14/40 → 233 kg/m³** (−35 %), **18/80 → 900 kg/m³** (×2.5) |
| Même table B' ? | **Oui, rigoureusement** — vérifié à `0.000e+00` par `bprime` |
| Même point de fonctionnement ? | **Non** — la pente `k = B'g/B'c` passe de 0.179 à 0.437 |
| ATG différentes ? | **Résine seule : identiques.** Composite : même forme, perte de 15.2 / 19.0 / 30.4 % |
| Paramètres d'Arrhenius ? | **Rien à extrapoler** — `A`, `E/R`, `m`, `f` inchangés ; seul le raccord au composite se ré-échelonne |
| Mise en données prête à l'emploi ? | § 5 — deux jeux complets (`F_i` ou `Δρ_i`), vérifiés équivalents |
| Le modèle est-il le bon ? | **Oui** — 18/43 et 18/33, réellement fabriqués, tombent exactement dessus (§ 6) |
| Réserve majeure | le rendement en char mesuré en jet plasma est **0.447**, pas 0.6198 (§ 6) |
| Porosité ? | **0.840 / 0.749 / 0.338** — c'est là que 18/80 cesse d'être un ablateur léger |

Le résultat structurant est une **séparation des rôles** :

- **XX ne fixe que l'échelle.** Il multiplie toutes les masses volumiques par
  un facteur, et rien d'autre.
- **YY porte toute la thermochimie de la réponse matériau** : couplage `k`,
  perte de masse en ATG, et — avec XX — la porosité.
- **Ni l'un ni l'autre ne touche la table B'.**

---

## 1. Densités et masses

| | 14/40 | 18/50 | 18/80 |
|---|---|---|---|
| masse fibres [kg/m³] | 140.0 | 180.0 | 180.0 |
| masse résine [kg/m³] | 93.3 | 180.0 | 720.0 |
| **ρ vierge [kg/m³]** | **233.3** | **360.0** | **900.0** |
| **ρ char [kg/m³]** | **197.8** | **291.6** | **626.2** |
| gaz libéré [kg/m³] | 35.5 | 68.4 | 273.8 |

Formules, avec `m_f = 10·XX` et `w = YY/100` :

```
m_r   = m_f · w/(1−w)          ρ_v = m_f + m_r
m_rc  = m_r · Y                ρ_c = m_f + m_rc          Y = 0.6198
```

**18/80 n'est plus un ablateur basse densité.** À 900 kg/m³ il pèse 2.1 fois
le ZURAM mesuré (419) et 3.2 fois le TACOT (280) — il quitte la classe des
ablateurs légers à préforme fibreuse. Ce n'est pas disqualifiant, mais cela
change la mission visée : un bouclier lourd, à faible recul, pour des flux très
élevés.

> **Incertitude de 15 % sur l'absolu.** Ces valeurs sont **nominales**. Le vrai
> ZURAM 18/50 pèse 419 kg/m³ et non les 360 que donne le calcul, parce que sa
> préforme pèse 205 kg/m³ et non les 180 nominaux (§ 3 de `resine_zuram.md`,
> anomalie ouverte n° 2 : *« we have never measured the weight and controlled
> the density of the CBCF 18-2000 preform »*). Tant que la préforme n'est pas
> pesée, compter ±15 % sur toute densité absolue. Les **rapports** entre
> variantes, eux, sont exacts.

---

## 2. La table B' est inchangée — vérifié

`bprime` ne consomme que trois compositions élémentaires, plus T, P et B'g
(`src/apps/bprime.cpp:304-330`) :

| Entrée | Origine | Dépend de XX/YY ? |
|---|---|---|
| bord de couche limite | air `O:0.21, N:0.79` | non |
| gaz de pyrolyse | **la résine seule** | non |
| char | fibres + résine carbonisée = carbone pur | non |

Le gaz est produit par la résine seule : changer XX/YY change la **quantité**
produite, pas la **composition** — et cette quantité est justement `B'g`, le
paramètre d'entrée balayé par la table. Le char reste `C:1.0` puisque fibres et
résine carbonisée sont toutes deux du carbone pur. Densité, porosité et
fractions volumiques sont des grandeurs **volumiques**, alors que le bilan de
surface est écrit sur des fractions massiques **normalisées**.

### Vérification numérique

Plutôt que de réutiliser le même XML — ce qui ne prouverait rien —
`zuram_variantes_bprime.py` **recalcule** pour chaque variante la composition
du gaz par fermeture `résine − char` depuis son propre bilan de masse :

| variante | résine [kg/m³] | char | gaz | C | H | O |
|---|---|---|---|---|---|---|
| 14/40 | 93.3 | 57.8 | 35.5 | 0.162679 | 0.725109 | 0.112212 |
| 18/50 | 180.0 | 111.6 | 68.4 | 0.162679 | 0.725109 | 0.112212 |
| 18/80 | 720.0 | 446.2 | 273.8 | 0.162679 | 0.725109 | 0.112212 |

Dispersion : **2.8e-17**. Les masses changent d'un facteur 8, les proportions
non. Les tables B' correspondantes, sur 284 points (1 atm, 500→4000 K,
B'g ∈ {0, 0.1, 0.5, 1.0}) :

| variante | écart max B'c | écart max h_w |
|---|---|---|
| 14/40 | **0.000e+00** | **0.000e+00** |
| 18/50 | **0.000e+00** | **0.000e+00** |
| 18/80 | **0.000e+00** | **0.000e+00** |

**Une seule table sert toutes les variantes.** C'est le même résultat que
`../tacot_bprime/autre_materiau.md` établissait pour TACOT vs CPh70.

---

## 3. Mais le point de fonctionnement change

En ablation stationnaire, le matériau impose sa propre droite :

```
B'g / B'c = ṁ_g/ṁ_c = (ρ_v − ρ_c)/ρ_c = k
```

et **k ne dépend que de YY**. En posant `w = YY/100` et `Y` le rendement en
char de la résine :

```
k = w(1−Y) / (1 − w(1−Y))
```

XX se simplifie — vérifié pour XX ∈ {10, 14, 18, 25, 40}, k identique à
1e-15 près.

| YY | 40 | 50 | 80 |
|---|---|---|---|
| **k** | **0.1794** | **0.2347** | **0.4372** |

Points fixes résolus sur la table, à 1 atm :

| T [K] | 14/40 B'c | B'g | 18/50 B'c | B'g | 18/80 B'c | B'g |
|---|---|---|---|---|---|---|
| 1500 | 0.1708 | 0.0306 | 0.1696 | 0.0398 | 0.1652 | 0.0722 |
| 2500 | 0.1783 | 0.0320 | 0.1782 | 0.0418 | 0.1773 | 0.0775 |
| 3000 | 0.1906 | 0.0342 | 0.1932 | 0.0454 | 0.2008 | 0.0878 |
| 3300 | 0.2128 | 0.0382 | 0.2185 | 0.0513 | 0.2381 | 0.1041 |
| 3600 | 0.3917 | 0.0703 | 0.4181 | 0.0981 | 0.5389 | 0.2356 |

Dans le plateau d'oxydation (T < 2500 K) les trois variantes sont à moins de
3 % l'une de l'autre : le soufflage pyrolytique y joue peu. L'écart se creuse
au genou de sublimation — à 3600 K, 18/80 est **38 % au-dessus** de 14/40.

### Le recul, lui, s'inverse

`ṡ = ṁ_c/ρ_c = B'c·ṁ_e/ρ_c` : à B'c comparable, c'est ρ_c qui décide.

| T [K] | 14/40 | 18/50 | 18/80 |
|---|---|---|---|
| 1500 | 1.00 | 0.67 | **0.31** |
| 3000 | 1.00 | 0.69 | 0.33 |
| 3600 | 1.00 | 0.72 | 0.43 |

**18/80 recule trois fois moins vite que 14/40** à flux identique, malgré un
B'c plus élevé. C'est le compromis classique : masse contre durée de vie.

---

## 4. Les ATG

### Sur la résine seule : identiques

Mêmes `A`, `E/R`, `m` ⇒ mêmes températures de pic, mêmes vitesses, même perte
finale de **38.02 %**. C'est le même polymère.

Les `f` du classeur (`'ZURAM_official'!A79:F83`) sont **déjà ramenées à la
résine seule** — la note `B77` précise que les `F` de la source, exprimées sur
le composite, ont été divisées par la fraction massique de résine. Elles sont
donc **transposables telles quelles** à toute variante.

| réaction | f | log₁₀A | E/R [K] | m |
|---|---|---|---|---|
| 1 | 0.035070 | 5.33 | 8 178.5 | 4.30 |
| 2 | 0.027687 | 8.69 | 16 068.4 | 3.70 |
| 3 | 0.095981 | 10.60 | 21 612.9 | 2.57 |
| 4 | 0.221495 | 11.67 | 26 423.8 | 4.63 |

### Sur le composite : même forme, autre amplitude

| variante | w_résine | perte totale | masse résiduelle |
|---|---|---|---|
| 14/40 | 0.400 | **15.21 %** | 84.79 % |
| 18/50 | 0.500 | **19.01 %** | 80.99 % |
| 18/80 | 0.800 | **30.42 %** | 69.58 % |

Les fibres Calcarb ne pyrolysent pas : elles **diluent** le signal. Une ATG de
composite ne distingue donc pas deux résines différentes d'une même résine en
proportion différente — il faut normaliser par `w_résine` avant toute
comparaison. Pour un code qui attend des `f` sur le composite, multiplier :

| variante | f₁ | f₂ | f₃ | f₄ |
|---|---|---|---|---|
| 14/40 | 0.014028 | 0.011075 | 0.038393 | 0.088598 |
| 18/50 | 0.017535 | 0.013843 | 0.047991 | 0.110748 |
| 18/80 | 0.028056 | 0.022150 | 0.076785 | 0.177196 |

### Les paramètres d'Arrhenius : rien à extrapoler

`A`, `E/R` et `m` décrivent la **chimie du novolac**, pas le composite. Ils
sont **identiques** pour toute variante. `f` aussi, dès lors qu'on l'exprime
sur la résine — ce que le classeur fait déjà. La seule opération est le
**raccord au composite**.

#### La forme normalisée, d'abord

Le classeur écrit ses paramètres dans la forme « TACOT-like » (note de version
1.4.0) :

```
d(ρ_i)/dt = − A_i · ρ_i,v · [(ρ_i − ρ_i,c)/ρ_i,v]^m · exp(−E_i/(R·T))
```

mais il ne fournit **que** `f_i, A_i, E_i/R, m_i` — jamais `ρ_i,v` ni `ρ_i,c`
séparément. Or cette équation n'est pas invariante d'échelle : en posant
`y = (ρ_i − ρ_i,c)/(ρ_i,v − ρ_i,c)`,

```
dy/dt = − A_i · [(ρ_i,v − ρ_i,c)/ρ_i,v]^(m−1) · y^m · exp(−E/RT)
```

Le facteur entre crochets **ne disparaît que si `ρ_i,c = 0`**. Une seule
lecture rend donc la donnée du classeur complète : chaque réaction consomme
intégralement sa pseudo-phase, le char étant porté par une phase inerte
distincte. Avec `x_i = ρ_i/ρ_i,v ∈ [0,1]` :

```
dx_i/dt = − A_i · x_i^m · exp(−E_i/(R·T))
ρ_résine(t) = ρ_résine,v · [ (1 − Σf_i) + Σ f_i·x_i ]
```

> Si votre code attend un `ρ_i,c` non nul, il faut corriger `A` du facteur
> `[(ρ_i,v − ρ_i,c)/ρ_i,v]^(m−1)`. Avec `m` entre 2.57 et 4.63, l'erreur est
> d'ordre un — ce n'est pas un détail.

#### (a) Base résine — inchangée pour toutes les variantes

| réaction | f [-] | log₁₀A | A [s⁻¹] | E/R [K] | m |
|---|---|---|---|---|---|
| 1 | 0.035070 | 5.33 | 2.138e+05 | 8 178.5 | 4.30 |
| 2 | 0.027687 | 8.69 | 4.898e+08 | 16 068.4 | 3.70 |
| 3 | 0.095981 | 10.60 | 3.981e+10 | 21 612.9 | 2.57 |
| 4 | 0.221495 | 11.67 | 4.677e+11 | 26 423.8 | 4.63 |

`Σf = 0.380234`, soit un rendement en char de `0.619766` — la valeur retrouvée
indépendamment par les fractions volumiques (§ 5 de `resine_zuram.md`).

#### (b) Base composite — seules les fractions se ré-échelonnent

`F_i = f_i · w` avec `w = YY/100` :

| variante | w | F₁ | F₂ | F₃ | F₄ | Σ |
|---|---|---|---|---|---|---|
| 14/40 | 0.400 | 0.014028 | 0.011075 | 0.038393 | 0.088598 | 0.152094 |
| 18/50 | 0.500 | 0.017535 | 0.013843 | 0.047991 | 0.110748 | 0.190117 |
| 18/80 | 0.800 | 0.028056 | 0.022150 | 0.076785 | 0.177196 | 0.304187 |

#### (c) Base volumique — pour un code façon `Pyrolysis model`

`ρ_i,v = f_i · ε_m · ρ_résine,intr` et `ρ_i,c = 0`, en kg/m³ de composite :

| variante | ε_m | ρ₁,v | ρ₂,v | ρ₃,v | ρ₄,v | char | fibres | Σ |
|---|---|---|---|---|---|---|---|---|
| 14/40 | 0.0710 | 3.27 | 2.58 | 8.96 | 20.67 | 57.84 | 140.00 | 233.33 |
| 18/50 | 0.1369 | 6.31 | 4.98 | 17.28 | 39.87 | 111.56 | 180.00 | 360.00 |
| 18/80 | 0.5475 | 25.25 | 19.93 | 69.11 | 159.48 | 446.23 | 180.00 | 900.00 |

La somme redonne exactement `ρ_vierge` de la § 1.

> **Ne pas mélanger les deux conventions de pyrolyse.** Le classeur suppose
> une densité intrinsèque de matrice *inchangée* avec perte de **volume**
> (anomalie ouverte n° 3), alors que la cinétique décrit une perte de
> **densité** à volume constant. Les deux donnent le même `ρ_char` composite,
> mais pas la même porosité intermédiaire.

#### Vérification

`zuram_cinetique.py` intègre **séparément** la base résine et, pour chaque
variante, la base composite avec ses propres `ρ_i,v` — deux calculs sans
référence l'un à l'autre. ATG simulée à 20 K/min :

| T [°C] | perte résine | 14/40 | 18/50 | 18/80 |
|---|---|---|---|---|
| 400 | 4.530 % | 1.812 % | 2.265 % | 3.624 % |
| 600 | 25.715 % | 10.286 % | 12.858 % | 20.572 % |
| 700 | 33.013 % | 13.205 % | 16.507 % | 26.411 % |
| 1000 | 37.239 % | 14.896 % | 18.620 % | 29.792 % |

Écart maximal entre les deux intégrations : **1.5e-14**. La relation
`(ρ_v − ρ)/ρ_v = w · perte_résine(T)` est donc exacte à tout instant, pas
seulement à l'asymptote.

Pic de DTG à **830.0 K (556.8 °C)**, identique pour toutes les variantes.

> **Queue algébrique — à connaître.** À 1400 K la perte simulée de la résine
> atteint 37.564 % et non les 38.023 % de `Σf`. Ce n'est pas une erreur
> d'intégration (résultat convergé en pas ; il ne dépend que de la température
> finale) mais une propriété du modèle : avec `m > 1`, `dx/dt = −A·x^m` donne
> `x ~ t^(−1/(m−1))`, qui ne s'annule jamais. Le rendement en char de 0.6198
> est une valeur **asymptotique** ; à 1400 K le modèle en rend 0.6244. À
> prendre en compte pour toute comparaison avec une ATG réelle, qui s'arrête
> à température finie.

### Réserve physique — c'est ici que le raisonnement peut casser

Tout ce qui précède suppose une résine **identiquement réticulée**. À 80 % en
masse, la résine remplit **62 % du volume poreux** : on n'a plus un film mince
autour des fibres (15 % à 18/50) mais des amas épais, dont la cuisson
— exothermie, évacuation de l'eau et de l'ammoniac issus de la HMTA — n'a
aucune raison d'être la même.

C'est précisément le point de Sykes sur le TACOT : le degré de cuisson gouverne
la quantité de composés piégés, donc la composition élémentaire, donc le
rendement en char. Un ZURAM 18/80 mal cuit à cœur aurait une **autre** résine
au sens de la table B'.

S'y ajoutent des effets de transport pendant la pyrolyse — craquage secondaire
des volatils sur un trajet plus long — que Torres-Herrador quantifie par les
nombres `Bi`, `Py_I` et `Py_II`.

> **14/40 est une extrapolation sûre** (résine plus diluée qu'à 18/50, film
> encore plus mince). **18/80 demande une ATG neuve** et, idéalement, une
> analyse élémentaire de la résine telle que cuite dans ce format.

---

## 5. Deux mises en données complètes, prêtes à recopier

Générées et vérifiées par `zuram_mise_en_donnees.py`. Les deux sont sur base
**composite** (fibres + résine) et décrivent le même matériau.

### Le bloc d'Arrhenius — commun aux deux jeux et aux trois variantes

```
dx_i/dt = − A_i · x_i^m_i · exp(−E_i/(R·T))          x_i(0) = 1
```

| i | A_i [s⁻¹] | E_i/R [K] | E_i [J/mol] | m_i |
|---|---|---|---|---|
| 1 | 2.137962e+05 | 8178.520143 | 68 000 | 4.30 |
| 2 | 4.897788e+08 | 16068.386634 | 133 600 | 3.70 |
| 3 | 3.981072e+10 | 21612.942201 | 179 700 | 2.57 |
| 4 | 4.677351e+11 | 26423.836403 | 219 700 | 4.63 |

La colonne `n` du classeur vaut 0 pour les quatre réactions : l'exposant qui
agit est `m`.

### JEU A — Arrhenius + F_i

```
ρ(T) = ρ_vierge · [ 1 − Σ_i F_i · (1 − x_i) ]
```

| | 14/40 | 18/50 | 18/80 |
|---|---|---|---|
| ρ_vierge [kg/m³] | 233.3333 | 360.0000 | 900.0000 |
| ρ_char [kg/m³] | 197.8449 | 291.5579 | 626.2318 |
| porosité vierge | 0.840265 | 0.749002 | 0.338378 |
| porosité char | 0.867251 | 0.801046 | 0.546556 |
| **F₁** | 0.014028037 | 0.017535047 | 0.028056075 |
| **F₂** | 0.011074766 | 0.013843458 | 0.022149533 |
| **F₃** | 0.038392523 | 0.047990654 | 0.076785047 |
| **F₄** | 0.088598131 | 0.110747664 | 0.177196262 |
| Σ F_i | 0.152093458 | 0.190116822 | 0.304186916 |

Contrôle : `1 − ΣF_i = ρ_char/ρ_vierge`, exact à 1e-15.

### JEU B — Arrhenius + Δρ_i

```
d(ρ_i)/dt = − A_i · Δρ_i · (ρ_i/Δρ_i)^m_i · exp(−E_i/(R·T))     ρ_i(0) = Δρ_i
ρ(T) = ρ_inerte + Σ_i ρ_i                       ρ_inerte = ρ_char
```

| | 14/40 | 18/50 | 18/80 |
|---|---|---|---|
| ρ_vierge [kg/m³] | 233.3333 | 360.0000 | 900.0000 |
| ρ_inerte [kg/m³] | 197.8449 | 291.5579 | 626.2318 |
| **Δρ₁** [kg/m³] | 3.273209 | 6.312617 | 25.250467 |
| **Δρ₂** | 2.584112 | 4.983645 | 19.934579 |
| **Δρ₃** | 8.958255 | 17.276636 | 69.106542 |
| **Δρ₄** | 20.672897 | 39.869159 | 159.476636 |
| Σ Δρ_i | 35.488474 | 68.442056 | 273.768224 |

Contrôle : `ρ_inerte + ΣΔρ_i = ρ_vierge`, exact.

> **Si le code demande un couple `(ρ_i,v ; ρ_i,c)`** : mettre `ρ_i,v = Δρ_i`
> **et `ρ_i,c = 0`**. Les `A` ci-dessus ne sont valides que pour `ρ_i,c = 0`.

### Équivalence des deux jeux — vérifiée

Intégration **indépendante** des deux systèmes, ATG à 20 K/min :

| T [°C] | 14/40 jeu A | jeu B | 18/50 jeu A | jeu B | 18/80 jeu A | jeu B |
|---|---|---|---|---|---|---|
| 400 | 229.1052 | 229.1052 | 351.8458 | 351.8458 | 867.3831 | 867.3831 |
| 600 | 209.3327 | 209.3327 | 313.7130 | 313.7130 | 714.8519 | 714.8519 |
| 900 | 199.0609 | 199.0609 | 293.9032 | 293.9032 | 635.6130 | 635.6130 |

Écart maximal sur toute la montée : **1.4e-11 kg/m³**.

Densités atteintes à 1400 K, à comparer au char asymptotique (queue algébrique,
cf. plus haut) : 198.27 / 292.38 / 629.54 contre 197.84 / 291.56 / 626.23.

> **Précision requise sur le rendement en char.** Le tronquer à
> `0.619766` au lieu de `0.619766355141076` suffit à désaccorder `ρ_char` et
> `ρ_vierge·(1 − ΣF)` dès le 7ᵉ chiffre. Les deux voies du classeur — fractions
> volumiques et somme des `f` — concordent en réalité à **3.3e-16**.

### Ce qui accompagne le bloc cinétique — identique aux trois

```
gaz de pyrolyse    C:0.171, H:0.722, O:0.107   (fractions molaires élémentaires)
char               C:1.0
bord de couche l.  O:0.21, N:0.79              (air)
table B'           la MÊME pour les trois       (vérifié à 0.000e+00, § 2)
```

Seules les densités, porosités et fractions volumiques distinguent les trois
mises en données. Tout ce qui est thermochimie de surface est commun.

## 6. Confrontation aux variantes réellement fabriquées

Pagan et al. (2017) publient cinq variantes de ZURAM testées en jet plasma à
l'IRS Stuttgart. Deux d'entre elles sont exactement l'exercice de ce document :
même préforme, même résine, teneur réduite.

| Variante | Préforme | Commentaire | ρ_pref | ρ_comp | ρ_char |
|---|---|---|---|---|---|
| 18/50 | carbone rigide | *standard quality* | 180 | 370 | 265 |
| 18/50-3 | carbone + CVI-SiC | revêtement SiC | 350 | 470 | 435 |
| **18/43** | carbone rigide | *resin content reduced by 25 %* | 180 | 340 | 250 |
| **18/33** | carbone rigide | *resin content reduced by 50 %* | 180 | 280 | 225 |
| M/50 | mullite | mullite au lieu de carbone | 128 | 360 | 230 |

### Ce que ça valide

La lecture de YY comme teneur en résine en % masse est **confirmée
exactement** : −25 % de résine donne 42.857 % → `18/43`, et −50 % donne
33.333 % → `18/33` (cf. `resine_zuram.md` § 1). Le modèle de ce document est
donc bien celui du fabricant.

### Ce que ça recadre

**Le DLR n'a exploré que la réduction de résine, jamais l'augmentation.**
Plage réelle : YY de 33 à 50.

- **14/40** — YY = 40 tombe *dans* la plage explorée : l'extrapolation est
  soutenue par deux matériaux réels. En revanche XX = 14 suppose une préforme
  Calcarb CBCF 14, absente de la famille — toutes gardent la CBCF 18.
- **18/80** — YY = 80 est *hors* de la plage, et dans la direction que le DLR
  n'a jamais prise. La réserve du § 4 cesse d'être une précaution théorique.

### Ce que ça contredit — le rendement en char

C'est le point le plus important de cette confrontation. Les densités de
Pagan impliquent, pour le 18/50 :

```
résine vierge     370 − 180 = 190 kg/m³
résine charbonnée 265 − 180 =  85 kg/m³
rendement en char  85/190   = 0.4474
```

La légende de la table confirme que c'est bien la logique employée
(*« scaled for remaining variants assuming identical ratios of solid resin
residue »*) : avec 0.4474 on retrouve 251.6 pour 18/43 (table : 250) et 224.7
pour 18/33 (table : 225).

**Or le classeur VKI donne 0.6198.** Avec cette valeur, ρ_char(18/50) vaudrait
`180 + 190·0.6198 = 297.8` kg/m³ au lieu des 265 mesurés — **12 % d'écart**.

L'explication la plus probable, à confirmer : les deux grandeurs ne décrivent
pas le même état. Le 0.6198 vient d'une ATG arrêtée vers 800–1000 °C ; le char
de Pagan est celui d'échantillons passés en jet plasma, bien au-delà de 2000 K.
Le modèle cinétique lui-même garde une queue algébrique au-delà de l'ATG
(§ 4) : un char plus chaud est plus léger.

**Conséquence pratique.** Pour une réponse matériau en ablation, le 0.4474 est
sans doute plus représentatif. Le couplage en dépend directement :

| rendement en char | origine | k = B'g/B'c |
|---|---|---|
| 0.6198 | classeur VKI (ATG) | 0.2410 |
| **0.4474** | Pagan, 18/50 mesuré | **0.3962** |

La table B', elle, reste inchangée : elle ne consomme pas le rendement en char
(§ 2). Mais toutes les densités et tous les `Δρ_i` du § 5 sont bâtis sur
0.6198 — à rebâtir sur 0.4474 si l'on privilégie la mesure post-plasma.

Vérification : `python zuram_bprime/nomenclature_pagan.py`

---

## 7. La porosité

| | 14/40 | 18/50 | 18/80 |
|---|---|---|---|
| ε fibres | 0.0888 | 0.1141 | 0.1141 |
| ε résine | 0.0710 | 0.1369 | 0.5475 |
| **porosité vierge** | **0.8403** | **0.7490** | **0.3384** |
| **porosité char** | **0.8673** | **0.8010** | **0.5466** |
| pores de la préforme remplis | 7.8 % | 15.5 % | **61.8 %** |

Trois conséquences pour 18/80 :

1. **Limite de fabrication.** Pores pleins, XX = 18 plafonne à **YY = 86.6 %**.
   80 est donc atteignable en principe, mais à 62 % de remplissage on est loin
   de l'infiltration légère du 18/50 — il faudrait vérifier que le moule RTM
   du DLR y arrive sans porosité fermée ni gradient.
2. **Perméabilité effondrée.** À porosité 0.34 contre 0.75, l'évacuation des
   gaz de pyrolyse devient un problème : surpression interne, délaminage.
   Le classeur donne 1.06e-11 m² pour le ZURAM vierge. Une loi de
   Carman-Kozeny en `ε³/(1−ε)²` — indication d'ordre de grandeur, pas une
   mesure — donne un facteur **75, soit 1.9 décade** entre 18/50 et 18/80
   (6.67 → 0.0885), contre un facteur 3.5 en sens inverse pour 14/40.
3. **Conductivité en hausse.** Moins de pores, plus de matière : les lois
   ajustées du classeur (`D43:E47`) ne sont plus applicables telles quelles.

À l'inverse, 14/40 reste dans le domaine où le ZURAM a été caractérisé, avec
une porosité un peu supérieure — extrapolation modérée.

---

## 8. Ce qu'il faudrait mesurer avant d'y croire

| Variante | Table B' | Réponse matériau | À mesurer |
|---|---|---|---|
| **14/40** | réutilisable telle quelle | fiable | densité de la préforme CBCF 14 |
| **18/80** | réutilisable **si** la résine cuit pareil | à valider | ATG + analyse élémentaire de la résine ; perméabilité ; conductivité |

Dans les deux cas, la première mesure à faire est celle qui manque déjà pour
le 18/50 : **peser la préforme**. C'est l'anomalie ouverte n° 2, et elle porte
14 % d'erreur sur toutes les densités absolues.

---

## 9. Fichiers

| Fichier | Rôle |
|---|---|
| `zuram_variantes.py` | bilan de phase solide, ATG, séparation des rôles |
| `zuram_variantes_bprime.py` | vérification par `bprime` : table identique, points de fonctionnement |
| `zuram_cinetique.py` | transposition des paramètres d'Arrhenius, ATG simulée |
| `zuram_mise_en_donnees.py` | les deux jeux prêts à recopier, et leur équivalence |
| `nomenclature_pagan.py` | confrontation à Pagan et al. 2017 : nom, densités, rendement en char |
| `resine_zuram.md` | traçabilité de la résine et de la nomenclature |
| `../tacot_bprime/autre_materiau.md` | même démonstration pour TACOT vs CPh70 |

Les deux scripts acceptent des variantes en argument :

```bash
python zuram_variantes.py 12/35 18/50 22/70
python zuram_variantes_bprime.py 14/40 18/80
```
