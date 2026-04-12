# Table B' du carbone graphite : démarche physique

## 1. Contexte — ablation d'un matériau carboné

Lors d'une rentrée atmosphérique, le bouclier thermique est soumis à un flux
de chaleur intense. Les matériaux carbonés (graphite, PICA, TACOT…) se
dégradent par **ablation** : le solide réagit chimiquement avec le gaz chaud
de la couche limite et/ou se sublime, injectant de la masse dans l'écoulement.

La quantification de ce débit de masse ablaté est l'objet de la **table B'**
(prononcer « B-prime »), outil standard en aérothermodynamique d'entrée
atmosphérique.

---

## 2. Définition des paramètres B'

La théorie de la couche limite laminaire conduit à deux flux de masse
adimensionnés définis par rapport au flux de masse du bord de couche limite
$\dot{m}_e$ :

$$
B'_c = \frac{\dot{m}_c}{\dot{m}_e}
\qquad
B'_g = \frac{\dot{m}_g}{\dot{m}_e}
$$

- $B'_c$ : débit de masse du **char** (fraction carbonée ablatée)
- $B'_g$ : débit de masse du **gaz de pyrolyse** (gaz issu de la décomposition
  interne du matériau)
- $\dot{m}_e$ : flux de masse au bord de la couche limite (référence)

Pour du graphite pur sans pyrolyse interne, $B'_g = 0$ et seul $B'_c$ est
non nul.

---

## 3. Bilan de masse élémentaire à la paroi

En supposant la couche limite en régime laminaire à nombre de Lewis unité, le
bilan de conservation de chaque **élément chimique** $i$ à la paroi s'écrit :

$$
(1 + B'_c + B'_g)\, Y_{w,i}
= Y_{e,i} + B'_c\, Y_{c,i} + B'_g\, Y_{g,i}
$$

où :

| Symbole | Signification |
|---------|---------------|
| $Y_{w,i}$ | Fraction massique de l'élément $i$ dans le gaz à la paroi |
| $Y_{e,i}$ | Fraction massique de l'élément $i$ au bord de couche limite (air) |
| $Y_{c,i}$ | Fraction massique de l'élément $i$ dans le char (graphite : $Y_{c,C}=1$) |
| $Y_{g,i}$ | Fraction massique de l'élément $i$ dans le gaz de pyrolyse |

En isolant $B'_c$ pour l'élément carbone ($i = C$) :

$$
\boxed{
B'_c = \frac{Y_{e,C} + B'_g\, Y_{g,C} - Y_{w,C}(1 + B'_g)}
             {Y_{w,C} - Y_{c,C}}
}
$$

Pour du graphite pur ($Y_{c,C} = 1$) dans l'air ($Y_{e,C} = 0$) sans pyrolyse
($B'_g = 0$) :

$$
B'_c = \frac{-Y_{w,C}}{Y_{w,C} - 1} = \frac{Y_{w,C}}{1 - Y_{w,C}}
$$

La seule inconnue est donc $Y_{w,C}$, la fraction massique de carbone dans le
**gaz d'équilibre à la paroi**.

---

## 4. Hypothèse d'équilibre thermochimique à la paroi

L'hypothèse centrale de la table B' est que la composition du gaz à la paroi
est **à l'équilibre thermochimique** à la température $T_w$ et à la pression
$P$ :

$$
\text{composition à la paroi} = f(T_w,\, P,\, \text{composition élémentaire à la paroi})
$$

La composition élémentaire à la paroi est dominée par le char : on simule un
**char infini** en ajoutant une grande quantité de carbone élémentaire à la
composition du bord de couche limite, puis on normalise. L'équilibre est alors
calculé sur ce mélange carbone-dominé.

### Algorithme numérique (implémenté dans MutationPP)

1. Initialiser la composition élémentaire à la paroi :

   $$
   X_{w,i} = Y_{e,i} + B'_g\, Y_{g,i}
   $$

2. Ajouter une grande quantité $\delta = \max(200,\; 100\, B'_g)$ de carbone :

   $$
   X_{w,C} \mathrel{+}= \delta, \qquad \text{puis normaliser}
   $$

3. Convertir les fractions massiques élémentaires en fractions molaires
   élémentaires.

4. Calculer la **composition d'équilibre en phase gazeuse** à $(T_w, P)$ via
   le solveur multi-phases de MutationPP
   (`equilibriumComposition(..., IN_PHASE)`).

5. Calculer $Y_{w,C}$ à partir des fractions molaires des espèces gazeuses :

   $$
   Y_{w,C} = \frac{\sum_j \nu_{Cj}\, X_{w,j}\, M_C}
                  {\sum_j X_{w,j}\, M_j}
   $$

   où $\nu_{Cj}$ est le nombre d'atomes de C dans l'espèce $j$ et $M_j$ sa
   masse molaire.

6. En déduire $B'_c$ par la formule du §3.

7. Calculer l'enthalpie de paroi :

   $$
   h_w = \frac{\sum_j X_{w,j}\, h_j(T_w)}{\sum_j X_{w,j}\, M_j}
   $$

---

## 5. Rôle de la phase condensée : C(gr)

La phase graphite condensée **C(gr)** doit être incluse dans le mélange. Sans
elle, le solveur d'équilibre traite tout le carbone comme gaz à toutes
températures : à 300 K, $\text{C}_3$ gazeux domine et $Y_{w,C} \approx 1$,
ce qui donne $B'_c \approx 200$ — valeur non physique.

Avec C(gr) dans le mélange :

- En dessous de ~3 900 K (à 1 atm), le graphite reste solide ; seule une
  faible fraction gazeuse (CO, CO₂…) est à l'équilibre → $Y_{w,C} \ll 1$
  → $B'_c$ faible et physiquement cohérent.
- Au-dessus de ~3 900 K, le graphite se sublime entièrement ; la phase gaz
  est dominée par C, C₂, C₃ → $B'_c$ élevé.

Le solveur multi-phases de MutationPP (`MultiPhaseEquilSolver`) gère
automatiquement la coexistence des phases gaz et solide.

---

## 6. Dépendance en pression

La table B' dépend de la pression $P$ car :

- La composition d'équilibre du gaz à la paroi est fonction de $(T_w, P)$.
- La température de sublimation du graphite augmente avec la pression
  (effet Clausius-Clapeyron).

À basse pression (0,001 atm), le graphite se sublime à plus basse température
et $B'_c$ monte plus tôt. À haute pression (1 000 atm), la sublimation est
retardée mais les réactions d'oxydation (C + O → CO) sont favorisées.

La table est donc une **surface 2D** : $B'_c = f(T_w,\, P)$, tracée comme
un faisceau de courbes isobares.

---

## 7. Mélange utilisé : `carbon-air.xml`

| Espèces gazeuses | Rôle |
|-----------------|------|
| N, O, NO, N₂, O₂ | Air de la couche limite |
| C, C₂, C₃ | Formes gazeuses du carbone |
| CO, CO₂ | Produits d'oxydation du carbone |
| CN | Produit de réaction C + N |

| Espèce condensée | Rôle |
|-----------------|------|
| C(gr) | Graphite solide — phase stable à basse T |

**Compositions nommées** (dans le fichier XML) :

- `air` : $N : 0{,}79,\ O : 0{,}21$ (fractions molaires élémentaires)
  → bord de couche limite
- `pyro` : $C : 1{,}0$ → gaz de pyrolyse pur carbone
  (utilisé pour $B'_g > 0$ et pour définir la composition du char)

---

## 8. Génération numérique de la table

Le binaire C++ `bprime` (MutationPP) est appelé en sous-processus par le
script Python `examples/python/carbon_bprime.py` pour chaque pression :

```
bprime -T 300:100:5000 -P <P_Pa> -b 0 -m carbon-air -bl air -py pyro
```

La sortie est une table ASCII : pour chaque $T_w$, la valeur de $B'_c$,
$h_w$ (MJ/kg), et les fractions molaires de toutes les espèces à la paroi.

Le script boucle sur **25 pressions** de $10^{-3}$ à $10^3$ atm
(espacement logarithmique) et :

1. Agrège les résultats dans `carbon_bprime_table.csv`.
2. Trace $B'_c(T_w)$ en échelle log₁₀ et $h_w(T_w)$ pour les 7 pressions
   puissances de 10, sauvegardé dans `carbon_bprime_table.png`.

---

## 9. Interprétation physique des résultats

| Régime | $T_w$ (1 atm) | Mécanisme dominant |
|--------|--------------|-------------------|
| Faible ablation | 300 – 700 K | Oxydation lente ; C(gr) stable |
| Oxydation active | 700 – 3 500 K | C(gr) + O → CO (réaction hétérogène) |
| Sublimation | 3 500 – 4 000 K | C(gr) → C, C₂, C₃ (gaz) |
| Sublimation totale | > 4 000 K | Plus de phase solide ; $B'_c \gg 1$ |

À haute pression, les mêmes régimes sont décalés vers des températures plus
élevées (la sublimation est retardée).

---

## 10. Réactions chimiques impliquées

Le calcul d'équilibre thermochimique à la paroi fait intervenir l'ensemble des
réactions entre les espèces du mélange. On les regroupe en trois familles.

### 10.1 Réactions hétérogènes (surface graphite ↔ gaz)

Ce sont les réactions entre le graphite solide C(gr) et les espèces gazeuses
de la couche limite. Elles sont à l'origine du débit d'ablation $B'_c$.

| # | Réaction | Régime dominant |
|---|----------|-----------------|
| R1 | $\text{C(gr)} + \text{O} \rightleftharpoons \text{CO}$ | 700 – 3 500 K (oxydation active) |
| R2 | $\text{C(gr)} + \tfrac{1}{2}\,\text{O}_2 \rightleftharpoons \text{CO}$ | 700 – 2 000 K |
| R3 | $\text{C(gr)} + \text{O}_2 \rightleftharpoons \text{CO}_2$ | 300 – 1 000 K (basse T) |
| R4 | $\text{C(gr)} + \text{CO}_2 \rightleftharpoons 2\,\text{CO}$ | > 1 000 K (Boudouard) |
| R5 | $\text{C(gr)} + \text{N} \rightleftharpoons \text{CN}$ | > 2 000 K |
| R6 | $\text{C(gr)} + \text{NO} \rightleftharpoons \text{CO} + \text{N}$ | 1 000 – 3 000 K |

> **Note :** dans l'approche équilibre thermochimique de la table B', ces
> réactions ne sont pas traitées individuellement avec des cinétiques ; elles
> contribuent collectivement à la composition d'équilibre du gaz à $(T_w, P)$.

### 10.2 Sublimation et dissociation du graphite

Au-delà de ~3 000 K, le graphite passe en phase gazeuse. Les formes gazeuses
stables dépendent de la température.

| # | Réaction | Régime |
|---|----------|--------|
| R7 | $\text{C(gr)} \rightleftharpoons \text{C}_{(g)}$ | > 3 500 K (atome C) |
| R8 | $2\,\text{C(gr)} \rightleftharpoons \text{C}_2{}_{(g)}$ | > 3 000 K |
| R9 | $3\,\text{C(gr)} \rightleftharpoons \text{C}_3{}_{(g)}$ | 2 500 – 4 000 K (molécule la plus stable à T intermédiaire) |

> **C₃** est la forme gazeuse du carbone thermodynamiquement la plus stable
> dans la plage 2 500 – 4 000 K. Sans la phase condensée C(gr), elle domine
> dès 300 K, ce qui donne des valeurs de $B'_c$ non physiques (cf. §5).

### 10.3 Réactions homogènes en phase gazeuse

Ces réactions se déroulent dans le gaz à la paroi et au sein de la couche
limite. Elles fixent l'équilibre entre les espèces C-N-O.

#### Dissociation / recombinaison de l'air

| # | Réaction |
|---|----------|
| R10 | $\text{N}_2 \rightleftharpoons 2\,\text{N}$ |
| R11 | $\text{O}_2 \rightleftharpoons 2\,\text{O}$ |
| R12 | $\text{N}_2 + \text{O}_2 \rightleftharpoons 2\,\text{NO}$ |
| R13 | $\text{N} + \text{O} \rightleftharpoons \text{NO}$ |

#### Chimie carbone–oxygène

| # | Réaction |
|---|----------|
| R14 | $\text{C} + \text{O} \rightleftharpoons \text{CO}$ |
| R15 | $\text{CO} + \text{O} \rightleftharpoons \text{CO}_2$ |
| R16 | $\text{CO}_2 \rightleftharpoons \text{CO} + \text{O}$ |
| R17 | $\text{C}_2 \rightleftharpoons 2\,\text{C}$ |
| R18 | $\text{C}_3 \rightleftharpoons \text{C}_2 + \text{C}$ |

#### Chimie carbone–azote

| # | Réaction |
|---|----------|
| R19 | $\text{C} + \text{N} \rightleftharpoons \text{CN}$ |
| R20 | $\text{CN} + \text{O} \rightleftharpoons \text{CO} + \text{N}$ |
| R21 | $\text{CN} + \text{O}_2 \rightleftharpoons \text{CO} + \text{NO}$ |

### 10.4 Synthèse par régime de température (1 atm)

```
T_w [K]    Phase stable    Espèce gazeuse dominante    Réactions clés
─────────────────────────────────────────────────────────────────────
 300 – 700    C(gr)          CO₂, O₂, N₂                R3
 700 –1500    C(gr)          CO, N₂                      R1, R2, R4
1500 –2500    C(gr)          CO, N, CN                   R1, R5, R10
2500 –3800    C(gr) + gaz    CO, CN, C₃                  R7–R9, R15
3800 –5000    Gaz seul       C, C₂, C₃                   R7–R9, R17–R18
```

---

## 11. Références

- Lees, L. (1956). *Laminar heat transfer over blunt-nosed bodies at hypersonic
  flight speeds*. Jet Propulsion.
- Kendall, R.M. et al. (1968). *An Analysis of the Coupled Chemically Reacting
  Boundary Layer and Charring Ablator*. NASA CR-1060.
- Scoggins, J.B. et al. (2020). *Mutation++: Multicomponent Thermodynamic and
  Transport Properties for Ionized Plasmas written in C++*. SoftwareX.
