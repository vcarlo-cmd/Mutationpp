# Classeur des mises en données — `mise_en_donnees_materiaux.xlsx`

Un onglet par matériau étudié. Chaque onglet reconstruit **étape par étape et en
formules vivantes** la chaîne qui mène de la donnée expérimentale publiée aux
trois nombres portés dans `data/mixtures/*.xml`.

```bash
python build_mise_en_donnees.py      # (re)génère le classeur
python verification_xlsx.py          # évalue toutes les formules et confronte aux XML
```

## Contenu

| Onglet | Matériau | Route de calcul du gaz de pyrolyse | Composition du XML |
|---|---|---|---|
| `Synthèse` | — | tableau récapitulatif, typologie des routes, recette | — |
| `Constantes` | — | masses molaires, air, motifs de résine, bilan de surface | — |
| `TACOT` | carbone/phénolique poreux | **A** — sommation atomique d'une spéciation moléculaire (Sykes 1967) | C:0.206, H:0.679, O:0.115 |
| `CPh70` | carbone/phénolique dense | **A**, héritée du TACOT (même résine) | C:0.206, H:0.679, O:0.115 |
| `ZURAM` | carbone/phénolique DLR | **B** — conversion masse → mole d'une analyse élémentaire mesurée | C:0.171, H:0.722, O:0.107 |
| `SC-1008` | résol (matrice du PICA) | **C** — fermeture élémentaire motif + rendement en char | C:0.2526, H:0.6407, O:0.1068 |
| `Liège-phénolique` | liège / phénolique 80/20 | **D** — fermeture sur DEUX constituants pyrolysants | C:0.287, H:0.592, O:0.121 |
| `Carbone` | graphite | **E** — stœchiométrie directe, pas de gaz | C:1.0 (char) |
| `Silice` | SiO₂ | **E** — char multi-élément, `-char-elem Si` | Si:1.0, O:2.0 |

## Trame commune à chaque onglet matériau

`§0` fiche d'identité · `§1` hypothèses · `§2` données d'entrée (avec sources) ·
`§3` calcul pas à pas · `§4` résultat porté dans le XML · `§5` contrôles et
validation croisée · `§6` sensibilité · `§7` réponse matériau (hors XML) ·
`§8` bloc XML.

## Conventions de couleur

| | |
|---|---|
| **bleu** | donnée d'entrée en dur (mesure, valeur publiée, hypothèse) — c'est ici qu'on édite |
| noir | cellule calculée par formule vivante |
| **vert** | lien vers un autre onglet |
| **rouge sur fond rosé** | la valeur telle qu'elle figure dans le fichier XML du dépôt |
| fond jaune | hypothèse clé |

Modifier une cellule bleue recalcule toute la chaîne — jusqu'à l'écart affiché
au `§4` avec la valeur réellement présente dans le dépôt.

## Vérification

`verification_xlsx.py` fait trois choses :

1. évalue **toutes** les formules du classeur (moteur `formulas`, pur Python) et
   signale toute erreur Excel ;
2. compare la composition calculée par chaque onglet à celle lue dans
   `data/mixtures/*.xml` — l'écart doit rester sous 0.1 %, soit l'arrondi des
   fichiers eux-mêmes ;
3. réinjecte les valeurs calculées dans le cache du fichier, pour qu'il s'ouvre
   déjà renseigné.

Sortie attendue :

```
543 formules évaluées, 0 erreur(s)
[OK ] TACOT              tacot_pyro       ... écart max 0.06 %
[OK ] CPh70              cph70_pyro       ... écart max 0.06 %
[OK ] ZURAM              VKIZuramPyroGas  ... écart max 0.03 %
[OK ] SC-1008            sc1008_pyro      ... écart max 0.01 %
[OK ] Liège-phénolique   cork_pyro        ... écart max 0.04 %
```

> LibreOffice n'est pas utilisé pour le recalcul : il est indisponible dans
> certains environnements de calcul. L'évaluation pure Python de
> `verification_xlsx.py` le remplace et couvre les mêmes formules.

## Sources

Ce classeur ne produit aucune donnée nouvelle : il met en formules ce que
documentent déjà, en prose, les fichiers suivants du dépôt.

| Onglet | Document de traçabilité |
|---|---|
| tous | `tacot_bprime/mise_en_donnees_xml.md` |
| TACOT | `tacot_bprime/resine_tacot.md` |
| CPh70 | `tacot_bprime/cph70_vs_tacot.md`, `cph70_bprime/README.md` |
| ZURAM | `zuram_bprime/resine_zuram.md` |
| SC-1008 | `sc1008_bprime/resine_sc1008.md` |
| Liège | `cork_bprime/mise_en_donnees_cork.md`, `cork_pyrolysis_data.py` |
| Carbone | `carbon_bprime/bprime_carbon_physique.md` |
