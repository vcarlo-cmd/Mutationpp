# Récapitulatif — données à renseigner dans le XML pour une table B'

Référence pratique pour créer `data/mixtures/<materiau>-air.xml`.
Tout ce qui suit est vérifié contre le parseur :
`src/general/MixtureOptions.cpp`, `src/thermo/Composition.cpp`,
`src/general/Mixture.cpp`.

---

## 1. Squelette

```xml
<mixture thermo_db="NASA-9">

    <species>
       <!-- liste des espèces GAZ + la ou les phases CONDENSÉES du char -->
       C H O N CH4 CN CO CO2 C2 C2H C2H2,acetylene C3 C4 C4H2,butadiyne C5
       HCN H2 H2O N2 CH2OH CNN CNC CNCOCN C6H6 HNC C(gr)
    </species>

    <element_compositions default="air">
        <composition name="air">N:0.79, O:0.21</composition>
        <composition name="xxx_pyro">C:0.206, H:0.679, O:0.115</composition>
        <composition name="xxx_char">C:1.0</composition>
    </element_compositions>

</mixture>
```

Commande associée :

```bash
bprime -T 300:25:4000 -P 101325 -b 0.5 \
       -m <materiau>-air -bl air -py xxx_pyro -char xxx_char -char-elem C
```

---

## 2. Les trois données à déterminer

`bprime` ne consomme que **trois compositions élémentaires** (plus T, P, B'g) —
cf. `src/apps/bprime.cpp:304-330`.

| # | Composition | Rôle | Option | Comment l'obtenir |
|---|-------------|------|--------|-------------------|
| 1 | bord de couche limite | `Y_e` | `-bl` | atmosphère : air `N:0.79, O:0.21` ; Mars `C:0.., O:..` ; etc. |
| 2 | gaz de pyrolyse | `Y_g` | `-py` | composition élémentaire du gaz produit par la **résine seule** |
| 3 | char | `Y_c` | `-char` + `-char-elem` | composition élémentaire du **résidu solide** |

> Sans `-char`/`-char-elem`, `bprime` utilise un char de carbone pur intégré.
> Les deux options vont **toujours ensemble** (erreur sinon).

### 2.1 Gaz de pyrolyse — à partir d'une composition moléculaire

Si vous disposez des fractions molaires des espèces du gaz, sommez les atomes :

```
n_E = Σ_j  ν_{E,j} · x_j          puis normaliser
```

Exemple TACOT (feuille `Pyrolysis-gas chemistry`) :

| Espèce | H2 | H2O | CH4 | C6H5OH | CO | CO2 | C6H6 |
|--------|-----|-----|-----|--------|-----|-----|------|
| x | 0.4992 | 0.2336 | 0.1000 | 0.0891 | 0.0576 | 0.0157 | 0.0048 |

→ C = 0.7368, H = 2.4291, O = 0.4117 mol → **C:0.206, H:0.679, O:0.115**

> Bloc de PYROLYSE (cinétique) du TACOT, dans ses trois écritures :
> `mise_en_donnees_pyrolyse.md`. Attention au facteur 9 sur la phase B.
>
> Traçabilité complète de ces sept nombres jusqu'à la table I de Sykes
> (NASA TN D-3810, 1967), et identification de la résine dont ils sortent :
> voir `resine_tacot.md` et `resine_tacot_verification.py`.
>
> Le ZURAM procède autrement : sa composition élémentaire est mesurée
> directement en fractions massiques, sans passer par une spéciation
> moléculaire. Voir `../zuram_bprime/resine_zuram.md`.

Ce qui donne, indifféremment (la normalisation est automatique, §3.1) :

```xml
<composition name="xxx_pyro">C:0.206, H:0.679, O:0.115</composition>
<composition name="xxx_pyro">C:0.7368, H:2.4291, O:0.4117</composition>
```

### 2.2 Char — cas mono-élément et multi-élément

**Char de carbone pur** (fibres de carbone + résine phénolique carbonisée) :

```xml
<composition name="xxx_char">C:1.0</composition>
```
avec `-char-elem C`.

**Char multi-élément** — si les constituants ne carbonisent pas vers le même
produit (fibres de silice, SiC, char retenant H/O…), il faut **pondérer par les
masses de char de chaque constituant** :

```
n_E = Σ_i  (m_char,i / M_i) · ν_{E,i}
```

Exemple silice pure (cf. `silice_bprime`) :

```xml
<composition name="silice">Si:1.0, O:2.0</composition>
```
avec `-char-elem Si`.

> `-char-elem` désigne l'élément sur lequel le bilan est résolu. Choisir un
> élément **présent dans le char et absent (ou minoritaire) au bord de couche
> limite** — C pour un char carboné dans l'air, Si pour la silice.

---

## 3. Règles du parseur — à connaître

### 3.1 Normalisation automatique

`Composition::componentsFromList` **normalise la somme à 1**. Vous pouvez donc
saisir indifféremment des moles brutes, des pourcentages ou des fractions.

### 3.2 `type="mole"` (défaut) ou `type="mass"`

```xml
<composition name="pyro" type="mass">C:0.51, H:0.09, O:0.40</composition>
```

Mutation++ convertit automatiquement (`XE_TO_YE` / `YE_TO_XE`,
`Mixture.cpp:131-138`). `bprime` demande des **fractions massiques** ; si vous
déclarez en mole, la conversion est faite pour vous.

> Piège : une composition en moles et la même en masse ne sont pas les mêmes
> nombres. Déclarez le `type` qui correspond à votre source.

### 3.3 Format et séparateurs

- format `nom:valeur, nom:valeur, ...`
- séparateurs acceptés : `:` `,` espace, tabulation, retour ligne
- nombre de jetons **pair** obligatoire, valeurs numériques
- noms **uniques** (doublon = erreur de parsing)
- tout élément cité doit exister dans le mélange, sinon `LogicError`

### 3.4 Attributs disponibles

| Balise | Attribut | Défaut | Rôle |
|--------|----------|--------|------|
| `<mixture>` | `thermo_db` | — | base thermo (`NASA-9`, `NASA-7`, `RRHO`) |
| | `mechanism` | — | mécanisme cinétique gaz |
| | `viscosity` | — | modèle de viscosité |
| | `thermal_conductivity` | — | modèle de conductivité |
| | `gsi_mechanism` | — | interaction gaz-surface |
| | `state_model` | — | modèle d'état (`ChemNonEq1T`…) |
| `<element_compositions>` | `default` | — | composition par défaut (doit exister) |
| `<composition>` | `name` | **requis** | nom appelé par `-bl` / `-py` / `-char` |
| | `type` | `mole` | `mole` ou `mass` |

Pour une table B', seul `thermo_db` est nécessaire.

---

## 4. Liste d'espèces

### 4.1 Inclure impérativement la phase condensée du char

`C(gr)` (ou `SiO2(L)`, etc.) **doit** figurer dans `<species>`. Sans elle, le
solveur multiphase ne peut pas condenser le carbone du char en excès :
Y_w,C → 1 et B'c sature sur une valeur non physique.

Ce n'est pas une espèce « gaz » de la table de référence — c'est une exigence du
solveur d'équilibre.

### 4.2 La liste conditionne le résultat

Écart mesuré contre la table de référence TACOT 3.0 :

| Liste | B'c moyen | médian | < 1 % |
|-------|-----------|--------|-------|
| 25 espèces de la référence | **0.433 %** | **0.011 %** | **89.2 %** |
| 35 espèces (liste étendue) | 2.409 % | 0.434 % | 61.3 % |

Pour **reproduire** une table publiée, reprendre sa liste d'espèces à
l'identique. Pour un calcul neuf, une liste plus riche est physiquement plus
complète — mais alors ne pas s'attendre à retrouver la table de référence.

### 4.3 Noms d'espèces

Ils doivent correspondre **exactement** à la base thermo, virgules comprises
(`C2H2,acetylene`, `C6H5OH,phenol`, `C4H2,butadiyne`). Vérifier :

```bash
grep -n "^HNC " data/thermo/nasa9.dat
build/src/apps/checkmix <mon-melange>     # liste espèces, éléments, Mw, phase
```

---

## 5. Ce qui n'entre PAS dans le XML

Le bilan de masse de surface est écrit sur des **fractions massiques
élémentaires**, donc normalisées. Les grandeurs **volumiques** du matériau n'y
ont aucune place :

| Donnée | Dans le XML ? | Où elle sert |
|--------|---------------|--------------|
| masse volumique ρ_v, ρ_c | **non** | solveur de réponse matériau ; récession `ṡ = B'c·ṁe/ρc` |
| porosité | **non** | conduction, perméabilité |
| ratio fibres/résine | **non**\* | fixe la **quantité** de gaz (donc B'g), pas sa composition |
| rendement en char | **non** | idem |
| cinétique de pyrolyse | **non** | fixe B'g au cours du temps |

\* sauf si les constituants donnent des chars de compositions élémentaires
**différentes** (fibres de silice + résine phénolique) : le ratio entre alors
dans `Y_c`, et la table change. Voir `autre_materiau.md`.

---

## 6. Recette pour un nouveau matériau

1. **Atmosphère** → composition `-bl` (air, Mars, Titan…).
2. **Résine** → composition moléculaire du gaz de pyrolyse → sommer les atomes →
   composition `-py`.
3. **Char** → si tous les constituants carbonisent vers le même élément :
   `C:1.0`. Sinon pondérer par les masses de char.
4. **Espèces** → tous les produits attendus des éléments présents,
   **plus la phase condensée**.
5. **Vérifier** : `checkmix` (noms, phases), puis un run court.
6. **Contrôler la physique** : à 300 K et B'g = 0 dans l'air, un char carboné
   doit donner B'c ≈ 0.087 (limite d'oxydation C + O₂ → CO₂, indépendante de la
   pression). Une valeur ≈ 200 signale l'absence de `C(gr)`.

---

## 7. Exemple complet commenté

```xml
<!-- Matériau X : composite carbone/phénolique, ablation dans l'air -->
<mixture thermo_db="NASA-9">

    <species>
       <!-- éléments et espèces C/H/O/N attendus à l'équilibre pariétal -->
       C H O N CH4 CN CO CO2 C2 C2H C2H2,acetylene C3 C4 C4H2,butadiyne C5
       HCN H2 H2O N2 CH2OH CNN CNC CNCOCN C6H6 HNC
       <!-- phase condensée : indispensable, sinon B'c non physique -->
       C(gr)
    </species>

    <element_compositions default="air">
        <!-- bord de couche limite (-bl) : air, fractions molaires -->
        <composition name="air">N:0.79, O:0.21</composition>

        <!-- gaz de pyrolyse (-py) : issu de la résine seule.
             Moles brutes acceptées, la normalisation est automatique. -->
        <composition name="x_pyro">C:0.206, H:0.679, O:0.115</composition>

        <!-- char (-char, avec -char-elem C) : fibres de carbone
             + résine carbonisée = carbone pur -->
        <composition name="x_char">C:1.0</composition>
    </element_compositions>

</mixture>
```

---

## 8. Récapitulatif en une ligne

> Dans le XML : **une liste d'espèces (avec la phase condensée)** et **trois
> compositions élémentaires** (bord de couche limite, gaz de pyrolyse, char).
> Rien d'autre. Toute l'information volumique du matériau vit dans le solveur de
> réponse matériau.
