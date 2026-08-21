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

## 5. La porosité

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

## 6. Ce qu'il faudrait mesurer avant d'y croire

| Variante | Table B' | Réponse matériau | À mesurer |
|---|---|---|---|
| **14/40** | réutilisable telle quelle | fiable | densité de la préforme CBCF 14 |
| **18/80** | réutilisable **si** la résine cuit pareil | à valider | ATG + analyse élémentaire de la résine ; perméabilité ; conductivité |

Dans les deux cas, la première mesure à faire est celle qui manque déjà pour
le 18/50 : **peser la préforme**. C'est l'anomalie ouverte n° 2, et elle porte
14 % d'erreur sur toutes les densités absolues.

---

## 7. Fichiers

| Fichier | Rôle |
|---|---|
| `zuram_variantes.py` | bilan de phase solide, ATG, séparation des rôles |
| `zuram_variantes_bprime.py` | vérification par `bprime` : table identique, points de fonctionnement |
| `resine_zuram.md` | traçabilité de la résine et de la nomenclature |
| `../tacot_bprime/autre_materiau.md` | même démonstration pour TACOT vs CPh70 |

Les deux scripts acceptent des variantes en argument :

```bash
python zuram_variantes.py 12/35 18/50 22/70
python zuram_variantes_bprime.py 14/40 18/80
```
