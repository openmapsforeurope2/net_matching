![under_construction](images/under_construction.png)



# Introduction

La présente documentation, à destination des développeurs, a pour objectif de présenter le détail du fonctionnement du processus de mise en cohérences des données de type réseau aux frontières ainsi que les principaux outils mis en oeuvre.

# Installation

## Code source 

Le code source de l'application est disponible sur le dépôt [area_matching](https://github.com/openmapsforeurope2/area_matching.git)

## Dépendances 

L'installation de l'application nécessite la compilation préalable de bibliothèques internes et externes à l'IGN.

Voici le graphe des dépendances :

<img src="images/dependencies.png" width="500" height="auto">

### Socle IGN 

Le socle logiciel de l'IGN regroupe un ensemble de bibliothèques développées en interne qui permettent d'unifier l'accès aux bibliothèques c++ de traitement et de stockage de données géographiques.
On y trouve notamment des modèles de données pivots (géométries, objet attributaire), des fonctions de lecture/écriture de conteneurs d'objets, des opérations sur les géométries, de nombreux algorithmes et outils spécifiquement conçus pour répondre à des problématiques géomaticiennes...

Le code source du socle ce trouve sur le dépôt [sd-socle](http://gitlab.forge-idi.ign.fr/socle/sd-socle.git)

### LibEPG 

Cette bibliothèque, développée à l'IGN et s'appuyant essentiellement sur le socle logiciel, contient de nombreux algorithmes et fonctions utilitaires dédiés spécifiquement aux besoins des produits européens (EGM/ERM) ainsi qu'au projet [OME2](https://github.com/openmapsforeurope2/OME2).
Elle comporte essentiellement des fonctions de généralisations, des fonctions utiles au management du processus tels que des utilitaires de log, d'orchestration, de gestion du contexte).
On y trouve également des opérateurs permettant d'encapsuler des objets géométriques complexes afin d'en optimiser la manipulation (par l'utilisation de graphes, d'indexes...) et ainsi d'accroitre les performances globales des processus.

Le code source de la bibliothèque libepg ce trouve sur le dépôt [libepg](http://gitlab.dockerforge.ign.fr/europe/libepg.git)


# Configuration

L'outil s'appuie sur de nombreux paramètres de configuration permettant d'adapter le comportement des algorithmes en fonctions des spécificités nationales (sémantique, précision, échelle, conventions de modélisation...).

On trouve dans le [dossier de configuration](https://github.com/openmapsforeurope2/net_matching/tree/main/config) les fichiers suivants :

- epg_parameters.ini : regroupe des paramètres de base issus de la bibliothèque libepg qui constitue le socle de développement l'outil. Ce fichier est aussi le fichier chapeau qui pointe vers les autres fichiers de configurations.
- db_conf.ini : informations de connexion à la base de données.
- theme_parameters.ini : configuration des paramètres spécifiques à l'application.


# Fonctionnement du processus

Le traitement de raccordement des objets linéaires est lancé pour un couple de pays frontaliers.

## Etapes préliminaires

Les données sur lesquelles ce traitement est lancé doivent avoir été nettoyées en amont à l'aide de l'outil **clean** du projet [data-tools](https://github.com/openmapsforeurope2/data-tools) qui permet de supprimer les objets trop éloignés de leur pays.
Cet outil doit être utilisé sur des tables de travail dans lesquelles sont extraites les données des deux pays à traiter autour de leur frontière commune.


## Principe général du traitement

Le processus de mise en cohérence des réseaux est décomposé en une succession d'étapes clés.
Afin d'orchestrer l'enchainement de ces étapes l'application utilise l'outil **epg::step::StepSuite** de la bibliothèque **libepg**. Ce dernier permet de lancer une succession de **epg::step::Step** dans lesquels sont implémentés les traitements de chaque étape.
Un code (numéro à trois chiffres) est attribué à chaque étape. Les étapes sont ordonnancées selon cette numérotation. Si une étape transforme les données sur lesquelles elle travaille, une ou plusieurs tables dédiées préfixées du code de l'étape sont créées. Ces créations sont réalisées en copiant les tables d'une étape antérieure (qui n'est pas nécessairement l'étape immédiatement antérieure, car toutes les étapes ne travaillent pas sur les mêmes données).
Ce fonctionnement permet de conserver les résultats intermédiaires du processus. Cela donne la possibilité d'arrêter et de reprendre le traitement en cours de processus et facilite le travail de d'analyse et de deboggage.


Les étapes qui composent le traitement de raccordement sont les suivantes :

**201** - initialisation du champ 'fictitious' des arcs du réseau
<br>
**202** - appariement des carrefours
<br>
**204** - généralisation des surfaces étroites en linéaires (les contours des surfaces sont constituées des arcs des 2 pays et les linéaires résultant de la fusion de ces arcs sont bi-nationaux)
<br>
**210** - génération des _'connectings lines'_ projetées aux frontières pour chacun des pays frontaliers
<br>
**211** - fusion des _'connecting lines'_ de chacun des pays en linéaires bi-nationaux
<br>
**212** - accrochage des _'connecting lines'_ dont les extrémités sont proches pour éviter les petites coupures
<br>
**213** - suppression des _'connecting lines'_ dont les couples de linéaires sont incohérents selon l'angle ou la distance
<br>
**214** - calcul de la géométrie des _'connecting lines'_ par interpolation
<br>
**220** - connection du réseau aux _'connecting lines'_
<br>
**230** - import des _'connecting lines'_ dans le réseau
<br>
**240** - génération des _'connecting points'_
<br>
**250** - connection du réseau aux _'connecting points'_
<br>
**255** - généralisation des surfaces étroites en linéaires bi-nationaux
<br>
**260** - nettoyage des artefacts (faces étroites, antennes, arcs superposés...)
<br>
**270** - connection des antennes hors pays au réseau du pays voisin
<br>
**280** - nettoyage des artefacts (faces étroites, antennes, arcs superposés, arcs de petite taille...)

> _Précisions_ :
> - _'Connecting line' : arc résultant de la fusion de deux arcs de deux pays différents réprésentant la même portion de réseau._
> - _'Connecting point' : sommet représentant un point de passage du réseau à la frontière. Les réseaux des deux pays limitrophes doivent être connectés à ce point afin d'assurer la continuité topologique du référentiel européen._
> - _champ 'fictitious' : la valeur de ce champ est 'vraie' lorsqu'un arc est couvert par une géométrie surfacique représentant le même objet (certains réseaux possèdes deux représentations linéaire et surfacique modélisées par deux classes d'objets différentes. Tous les objets linéaires n'ont pas de représentation surfacique)._

L'outil **epg::step::StepSuite** donne la possibilité de ne lancer que certaines étapes ou une plage de plusieurs étapes.

La liste de l'ensemble des étapes qui constituent le processus de raccordement diffère selon le réseau traité :
- cours d'eau (table hy.watercourse) : étapes 201 à 280
- tronçons routiers (table tn.road_link) :  étapes 210 à 280
- tronçons de voies ferrées (table tn.railway_link) : étapes 240 à 280

## Configuration

L'outil s'appuie sur de nombreux paramètres de configuration permettant d'adapter le comportement des algorithmes en fonctions des spécificités nationales (sémantique, précision, échelle, conventions de modélisation...).

On trouve dans le [dossier de configuration](https://github.com/openmapsforeurope2/net_matching/tree/main/config) les fichiers suivants :

- epg_parameters.ini : regroupe des paramètres de base issus de la bibliothèque libepg qui constitue le socle de développement l'outil. Ce fichier est aussi le fichier chapeau qui pointe vers les autres fichiers de configurations.
- db_conf.ini : informations de connexion à la base de données.
- theme_parameters_watercourse_link.ini : configuration des paramètres spécifiques au réseau de cours d'eau.
- theme_parameters_railway_link.ini : configuration des paramètres spécifiques au réseau de voies ferrées.
- theme_parameters_road_link.ini : configuration des paramètres spécifiques au réseau routier.


## Lancement du traitement

L'outil s'utilise en ligne de commande.
Le traitement peut être lancé sur trois types de réseaux :
- routier (code tn)
- férré (code ra)
- hydrographique (code hy)

<br>

Les paramètres sont les suivants :
* c [obligatoire] : chemin vers le fichier de configuration
* s [obligatoire] : suffix de la table de travail
* as [optionnel] : suffix des tables de travail des surfaces (à utiliser uniquement lors du traitement de la classe d'objets watercourse_link)
* t [obligatoire] : nom de la classe d'objet (doit être parmi les valeurs : road_link, railway_link, watercourse_link)
* sp [obligatoire] : code de l'étape(s) à executer (exemples: 220 ou 220,240 ou 210-280...)
* arguments libres [obligatoire] : codes des deux pays frontaliers

<br>

Exemple de lancement du traitement complet sur le couple de pays France (code pays 'fr') et Belgique (code pays 'be') pour le thème routier :
```
bin/net_matching --c path/to/config/epg_parameters.ini --s 20251118 --t road_link be fr
```

Exemple du lancement d'une seule étape :
```
bin/net_matching --c path/to/config/epg_parameters.ini --s 20251118 --t road_link -sp 240 be fr
```

Exemple de lancement d'une plage d'étapes :
```
bin/net_matching --c path/to/config/epg_parameters.ini --s 20251118 --t road_link -sp 240-260 be fr
```
Ici toutes les étapes de 240 à 260 (incluses) sont jouées.


## Les étapes - fonctionnement détaillé

### 201 : CorrectCountryConnectivity

Lors de cette étape il est question de corriger la cohérence topologique du réseau, car une bonne connectivité est nécessaire au bon déroulement du processus de raccordement. En effet, les données sources fournies par les organismes producteurs peuvent présenter des problèmes de connectivité (les arcs connectés doivent avoir un sommet en commun). Cela peut relever d'erreur lors de l'acquisition ou de choix délibérés de modelisation.
Par exemple, en ce qui concerne le réseau le réseau hydrographique, certain pays peuvent faire le choix de ne pas couper les tronçons des cours d'eaux principaux au niveau des points de confluence.

Deux types de corrections sont opérées par le processus:
- la création de connexion manquante par découpage des arcs situés à proximité de sommets d'autres arcs, puis connection du résultat de la découpe à ces sommets.
- la correction de la connexion de sommets très proches mais ne possédant pas précisemment les mêmes coordonnées.

#### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| EDGE_TABLE_INIT                | X      | X      | X                  | table du réseau à traiter                                   |

#### Principaux opérateurs de calcul utilisés :
- app::calcul::ConnectivityCorrectorOp

#### Description du traitement :
Paramètre utilisés: 
| paramètre                       | description                                                                                        |
|---------------------------------|----------------------------------------------------------------------------------------------------|
| CC_DIST_THRESHOLD               | distance seuil à partir de laquelle il faut établir une connexion                                  |
| NATIONAL_IDENTIFIER_NAME        | nom du champ pour l'identifiant national                                                           |


Pour la création des connexions manquantes, dans un soucis de performance, on charge préalablement en mémoire l'ensemble des géométries des arcs du réseau et on enregistre les relations d'adjacence existantes.
Ensuite on parcourt les sommets et pour chacun d'eux on detecte la présence éventuelle d'un arc non-adjacent 'e' situé à une distance inférieur à la distance seuil _CC_DIST_THRESHOLD_. Le cas échéant on enregistre ce point comme futur point de coupure de l'arc 'e'.
Enfin pour chacun des arcs possédant un ou plusieurs points de coupure, on procéde à la découpe. L'arc original est effacé et remplacé par les arcs issus de la découpe. De nouveaux identifiants _NATIONAL_IDENTIFIER_NAME_ sont créés pour les nouveaux arcs en prenant comme base l'identifiant de l'arc original auquel on ajoute un suffix incrémental (par exemple si un objet possédant l'identifiant <id> est découpé en trois objets ces objets auront pour identifiant <id>_1, <id>_2 et <id>_3).

![201_a](images/201_a_with_key.png)

Pour la corrections des connexion, comme précédemment, on charge l'ensemble des informations nécessaires au traitement en mémoire (géométries des arcs, informations sur les sommets...).
Le principe du traitement consiste à itérer sur les sommets et pour chaque sommet parcouru (qui devient le sommet de référence), chercher les autres sommets situés à une distance inférieures au seuil _CC_DIST_THRESHOLD_, puis, pour chacun des arcs auquels appartiennent ces sommets, on remplace sur les géométrie stockées en mémoire l'extrémité correspondante par la géométrie du point de référence.
L'ensemble des arcs modifiés sont enregistré en base de données en fin de traitement.
A noter que le point de référence (sommet restant fixe, vers lequel les autres sommets sont déplacés) est arbitrairement choisi parmi l'ensemble des sommets devant être connectés car cela dépend de l'ordre dans lequel l'algorithme parcourt les sommets.

![201_b](images/201_b_with_key.png)

### 202 : FillFictitiousField

Cette étape est spécifique au réseau hydrographique. Elle consiste à renseigner et/ou corriger le champ _EDGE_FICTITIOUS_NAME_ (indiquant si un arc est fictif) qui peut être manquant ou mal renseigné dans les donnnées sources livrées par les producteurs nationaux.
Le réseau hydrographique est réprésenté par deux classes d'objets linéaires (table watercourse_link) et surfaciques (tables watercourse_area et standing_water). Un arc du réseau d'objets linéaires est marqué comme _'fictif'_ s'il est recouvert par un objet surfacique représentant le même élément hydrographique du monde réel.

Il est à noter que la mise en cohérence des surfaces hydrographiques doit avoir été réalisée préalablement à l'initialisation du champ _EDGE_FICTITIOUS_NAME_ réalisé ici car le processus de mise en cohérence peut conduire à la suppression de surfaces et portions de surfaces.

#### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| EDGE_TABLE_INIT                | X      | X      | X                  | table du réseau à traiter                                   |
| WATERCOURSE_AREA_TABLE         | X      |        |                    | table des cours d'eau surfaciques                           |
| STANDING_WATER_TABLE           | X      |        |                    | table des surfaces d'eau                                    |


#### Principaux opérateurs de calcul utilisés :
- app::calcul::FillFictitiousFieldOp

#### Description du traitement :
Paramètre utilisés: 
| paramètre                       | description                                                                                        |
|---------------------------------|----------------------------------------------------------------------------------------------------|
| EDGE_FICTITIOUS_NAME            | nom du champ indiquant si l'objet est fictif                                                       |
| FFF_RATIO                       | ratio de recouvrement avec une surface à partir duquel un arc du réseau est considéré comme fictif |

On parcourt l'ensemble des arcs du réseau et pour chacun d'entre eux on calcule leur intersection avec les surfaces du même pays. On en déduit le ratio _longueur_recouverte/longeur_totale_. si ce ratio dépasse le seuil _FFF_RATIO_ le champ _EDGE_FICTITIOUS_NAME_ sera marqué comme 'vrai', sinon il sera marqué 'faux' (le champ peut avoir été injustement marqué 'vrai' dans les données sources).

![202](images/202.png)


### 203 : JunctionMatching

Dans cette étape il est question d'appairer les carrefours (noeuds de degré supérieur à deux) des deux pays. Puis, pour chaque couple de carrefours appairés, réaliser les déplacements nécessaires pour que leur localisation soit identique.
C'est une étape préparatoire qui permet de faciliter le travail de raccordement des zones complexes que sont les carrefours et d'éviter la création d'artefacts non désirés.


#### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| EDGE_TABLE_INIT                | X      | X      | X                  | table du réseau à traiter                                   |

#### Principaux opérateurs de calcul utilisés :
- app::calcul::JunctionMatchingOp

#### Description du traitement :
Paramètre utilisés: 
| paramètre                       | description                                          |
|---------------------------------|------------------------------------------------------|
| EDGE_FICTITIOUS_NAME            | nom du champ indiquant si l'objet est fictif         |
| JM_MAX_DIST                     | distance maximale entre deux carrefours appairables  |

Les réseaux des deux pays sont chargés dans deux graphs séparés.
La première étape de ce traitement consiste ensuite à parcourir tous les noeuds de degré supérieur à deux d'un graph et identifier leur meilleur candidat dans l'autre graph. L'opération inverse est ensuite réalisée. On peut alors dresser la liste des carrefours appairés (ce sont ceux qui sont réciproquement meilleurs candidats).
La seconde étape consiste à réaliser les déplacements pour que les deux carrefours convergent au même point. Si l'un des deux carrefours est fictif (possède au moins un arc incident fictif) le carrefour non-fictif est déplacé vers le carrefour fictif. Cela permet de conserver la cohérence entre les arc fictifs et les surfaces. Si les carrefours sont tous deux fictifs ou non-fictifs il sont tous les deux déplacés vers un point médian situé entre les deux carrefours.
Afin de minimiser les modifications géométriques sur les arcs du réseau seul l'extrémité des arcs sont déplacés, les points intermédiaires ne sont pas modifiés. Cela peut avoir pour conséquence la création de rebroussements.

![203_a](images/203_a.png)

Afin d'éviter la création de rebroussement le principe constiste à projeter le point cible (point vers lequel l'extrémité de la polyligne doit être déplacée) sur la polyligne à déplacer, réaliser une découpe de la polyligne avec le point projeté, puis déplacer l'extrémité de la polyligne découpé vers le point cible.

![203_b](images/203_b.png)
<br>
![203_c](images/203_c.png)


### 204 : GenerateCLinArea

L'objectif de cette étape et de fusionner les objets issus des modélisations des deux pays frontaliers représentant un même objet du monde réel. Les objets résultant de cette fusion sont des _connecting lines_. La valeur de ses attributs est calculée par concaténation des attributs des objets fusionnés.

#### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                                                                                  |
|--------------------------------|--------|--------|--------------------|------------------------------------------------------------------------------------------------------------------------------|
| EDGE_TABLE_INIT                | X      | X      | X                  | table du réseau à traiter (à renseigner si suffix (-s) non renseigné en ligne de commande)                                   |
| MATCHED_WATERCOURSE_AREA_TABLE | X      |        |                    | table des cours d'eau surfaciques raccordés aux frontières (à renseigner si suffix (-as) non renseigné en ligne de commande) |
| MATCHED_STANDING_WATER_TABLE   | X      |        |                    | table des surfaces d'eau  raccordés aux frontières (à renseigner si suffix (-as) non renseigné en ligne de commande)         |

#### Principaux opérateurs de calcul utilisés :
- app::calcul::CLInAreaGenerationOp

#### Description du traitement :
Paramètre utilisés: 
| paramètre                       | description                                                                                                                           |
|---------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| NATIONAL_IDENTIFIER_NAME        | nom du champ pour l'identifiant national                                                                                              |
| W_TAG_NAME                      | champ de travail permettant de marquer les objets traités                                                                             |
| CLA_SURFACE_WIDTH               | seuil de largeur des surfaces fines                                                                                                   |
| CLA_FICTITIOUS_RATIO_THRESHOLD  | ratio à partir duquel un chemin est considéré comme fictif (longueur des arcs fictifs/longeur totale)                                 |
| CLA_CL_LENGTH_THRESHOLD         | seuil de longeur en dessous duquel les connecting lines composant le contour d'une face fine peuvent être potentiellement effondrées  |
| CLA_CL_MIN_RATIO_IN_AREA        | ratio de longueur en dessous duquel les connecting lines composant le contour d'une face fine peuvent être potentiellement effondrées |
| LIST_ATTR_W                     | liste des attributs de travail (à ne pas fusionner)                                                                                   |
| LIST_ATTR_JSON                  | liste des attributs de type json (utilisé par l'opération de fusion des attributs)                                                    |


Le processus complet peut nécessiter de lancer plusieurs itérations de traitement. En effet, si à l'issu d'un traitement une ou plusieurs fusions ont été opérées alors une itération de traitement supplémentaire sera lancée, car la fusion peut générer de nouveaux cas à traiter. Le processus s'arrête lorsqu'à l'issu d'une itération aucune fusion n'a eu lieu.
Un traitement se déroule en plusieurs étapes:
- ré-initialisation du champ _W_TAG_NAME_ : ce champ est utilisé de manière interne au traitement afin de permettre aux opérateurs de savoir quels objets ont été antérieurement traités par d'autres opérateurs
- création d'un graph planaire à partir de l'ensemble des réseaux des deux pays frontaliers : cette étape permet le chargement en mémoire des données de travail, la création des faces, facilite la détection des arcs superposés...
- fusion des arcs des deux pays superposés
- fusion des arcs des deux pays constituant des faces fines : les faces fines sont traitées si celles-ci sont constituées de deux chemins appartenant aux deux pays et que les éventuelles connecting lines constituants ces chemins ont pu être substituées par un arc non fusionné (cette substitution se base sur le _NATIONAL_IDENTIFIER_NAME_)
- modification et objets possédant des portions fusionnée : seule une ou plusieurs portion(s) de la géométrie d'un arc peu(ven)t avoir été fusionnées, l'arc peut alors être découpé en plusieurs objets fusionnés et non-fusionnés
- modification (déplacement) des arcs incidents : les arcs incidents aux arcs fusionnés deviennent incidents aux _connecting lines_ résultant de la fusion, ils sont positionnés par conservation de l'abscisse curviligne. seule l'extrémité de l'arc incident est modifiée (pas de déformation amortie de l'arc incident)
- nettoyage des connecting lines superposées : correction des artefacts de traitement
- concaténation des arcs possédant le même _NATIONAL_IDENTIFIER_NAME_ : cette étape permet de corriger les découpages superflus créés par le traitement
- effondrement (suppression) des connecting lines : pour chaque face fine on effondre les connecting lines si celles-ci réprésentent une portion de son contour inférieure à _CLA_FICTITIOUS_RATIO_THRESHOLD_ et une longeur inférieure à _CLA_FICTITIOUS_LENGTH_THRESHOLD_. En cas d'effondrement, les éventuels arcs incidents sont déplacés avec un opérateur de déformation amortie.

A noter que la création des connecting lines dans les faces fines prend en compte la propriété _EDGE_FICTITIOUS_NAME_ afin de déterminer leur géométrie. Si l'un des deux objets fusionnés est fictif et l'autre non, c'est la géométrie de l'objet fictif qui sera prise, si les deux objets sont tout deux fictifs ou non-fictifs une géométrie moyenne sera calculée. En cas de fusion de deux arcs fictifs, on vérifie que le résultat de la fusion est bien entièrement inclus dans une surface fusionnée (résultant du traitement _net_area_matching_), si ce n'est pas le cas cette fusion est abandonnée.

![204_a](images/204_a_with_key.png)
<br>
![204_b](images/204_b_with_key.png)


### 210-214 : Generate connecting lines on Border

L'objectif de cette suite d'étapes et la génération de _connecting lines_ résultant de la fusion d'arcs des deux pays localisés à proximité des frontières.
Les _connecting lines_ sont enregistrées dans une table dédiée. C'est lors d'une étape postérieure que ces objets seront injectés dans la table des arcs en replacement des arcs dont ils sont le résultat de la fusion.
_
#### 210 : GenerateConnectingLinesByCountry

L'objectif de cette étape et de générer, pour chaque pays, les connecting lines correspondant à des arcs ou portions d'arcs longeant la frontière. 

##### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| CL_TABLE                       |        | X      | X                  | table des connecting lines                                  |
| EDGE_TABLE_INIT                | X      |        |                    | table du réseau à traiter                                   |

##### Principaux opérateurs de calcul utilisés :
- app::calcul::CFeatGenerationOp

##### Description du traitement :

Paramètre utilisés: 
| paramètre                     | description                                                                                                                           |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| CL_BUFFER_DIST                | rayon du buffer autour des frontières                                                                                                 |
| CL_THRESHOLD_NO_CL            | seuil de longueur de portion d'arc hors buffer determinant une interruption de connecting line                                        |
| CL_RATIO_IN_BUFFER            | proportion d'une portion d'arc qui doit être situé dans le buffer de la frontière pour pouvoir donner naissance à une connecting line |
| CL_SNAP_ON_VERTEX_BORDER_DIST | distance de snapping de la projection des extremités des arcs sur les points des polylignes des frontières                            |
| CL_BORDER_MAX_ANGLE           | angle maximum entre la frontièere et un arc ou une portion d'arc pour qu'ils puissent générer une connecting line                     |
| CFG_BOUNDARY_ANGLE_THRESHOLD  | seuil d'angle pour le découpage des frontières au niveau des angles aigu (si le paramètre vaut 180, aucun découpage n'est réalisé)    |
| LINKED_FEATURE_ID             | nom du champ indiquant l'identifiant de l'arc associé à une connecting line                                                           |

L'objectif, dans un premier temps, est de déterminer les arcs ou portions d'arcs qui donneront naissance à des connecting lines,
Pour cela un buffer est créé autour de la frontière et chaque arc intersectant celui-ci est découpé selon le contour du buffer.
Pour chaque arc on parcourt l'ensemble de sous-arcs ainsi obtenu et on détermine quelles sont les ensembles de sous-arcs adjacents pouvant donner naissance à des connecting lines.
Un ensemble de sous-arcs est éligible à la génération d'une connecting line si celui-ci n'est composé que de parties localisées à l'interieur du buffer et éventuellemnent de partie hors buffer de longueur inférieur au seuil _CL_THRESHOLD_NO_CL_. Un ensemble éligible donnera effectivement naissance à une connecting line si la proportion (en longueur) des parties situées dans le buffer est supérieur au seuil _CL_RATIO_IN_BUFFER_.

![210](images/210_with_key.png)


#### 211 : MergeConnectingLinesOnBorder

Le but est ici de fusionner deux à deux les connecting lines ou portions de connecting lines dont les arcs associés peuvent être appairés. A l'issu de ce traitement il ne doit subsister que des connectings associés à deux arcs des deux pays.

##### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| CL_TABLE                       |        | X      | X                  | table des connecting lines                                  |

##### Principaux opérateurs de calcul utilisés :
- app::calcul::CFeatGenerationOp

##### Description du traitement :

Paramètre utilisés: 
| paramètre                   | description                                                                                                               |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------|
| CL_EDGE_MAX_DIST            | distance de hausdorff maximum entre deux portions d'arc pour quelles puissent être fusionnées en une connecting line      |
| CL_SNAP_PROJ_CL_2_EDGE_DIST | distance de snapping de la projection des extremités des connecting lines sur les points des polylignes des arcs associés |
| LINKED_FEATURE_ID           | nom du champ indiquant les identifiants des arcs associés à une connecting line                                           |


La première étape consiste à construire un graph planaire à partir des connecting lines. Les linéaires superposés sont ainsi agrégés et découpés pour former le graph, chaque connexion du graph est ainsi associée à un ou plusieurs arcs du réseau traité.
Pour les connexions du graph possédant plusieurs arcs associés provenant des deux pays, on cherche les meilleurs couples d'arcs à appairer.
Pour qu'un couple d'arcs soit appairable il faut que :
- les deux arcs soient de pays différents
- la distance de hausdorff entre les portions d'arcs associés à la connecting line soit inférieure au seuil _CL_EDGE_MAX_DIST_

Pour chaque couple appairé une connecting line est créée. Toutes les connecting lines ou portions de connecting line non fusionnées sont supprimées.

![211_a](images/211_a_with_key.png)
<br>
![211_b](images/211_b_with_key.png)


#### 212 : SnapConnectingLines

Cette étape vise à corriger les discontinuités entre les connecting lines qui ont pu être crées par les traitements précédents.

##### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| CL_TABLE                       |        | X      | X                  | table des connecting lines                                  |

##### Principaux opérateurs de calcul utilisés :
- app::calcul::CFeatGenerationOp

##### Description du traitement :

Paramètre utilisés: 
| paramètre                   | description                                                                                                      |
|-----------------------------|------------------------------------------------------------------------------------------------------------------|
| CL_CL_CLOSEST_MAX_DIST      | seuil de distance entre les extrémités de deux connecting lines en dessous duquel une connection doit être créée |

Pour détecter les discontinuités, un graph des connecting lines est chargé. Pour chaque noeud de degré 1 on recherche s'il existe des extrémités d'autres connecting lines situés à une distance inférieure au seuil _CL_CL_CLOSEST_MAX_DIST_. 
Les extremités proches sont déplacées vers le barycentre de l'ensemble de ces points.

![212](images/212_with_key.png)



#### 213 : DeleteConnectingLines

Cette étape consiste à supprimer les connecting lines pour lesquelles les arcs associés ne peuvent être fusionnés car trop éloignés ou parce que présentant des orientations trop différentes.

##### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| CL_TABLE                       |        | X      | X                  | table des connecting lines                                  |
| EDGE_TABLE_INIT                | X      |        |                    | table du réseau à traiter                                   |

##### Principaux opérateurs de calcul utilisés :
- app::calcul::CFeatGenerationOp

##### Description du traitement :

Paramètre utilisés: 
| paramètre                   | description                                                                                                               |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------|
| CL_EDGE_MAX_ANGLE           | angle maximum autorisé entre les deux portions d'arcs associées à une connecting line                                     |
| CL_EDGE_MAX_DIST            | distance de hausdorff maximum autorisée entre les deux portions d'arc associées à une connecting line                     |
| CL_MIN_LENGTH               | longueur minimum d'une connecting line isolée (non connectée à d'autres connecting lines)                                 |
| CL_SNAP_PROJ_CL_2_EDGE_DIST | distance de snapping de la projection des extremités des connecting lines sur les points des polylignes des arcs associés |
| LINKED_FEATURE_ID           | nom du champ indiquant les identifiants des arcs associés à une connecting line                                           |

On charge dans un premier un graphe des connecting lines, puis, pour chaque connecting line dont au moins une des extrémités n'est pas connectée à d'autres connecting lines (noeud de degré 1), on calcule la portion pour chacun des deux arcs associés correspondant à la projection de la connecting line sur l'arc. On obtient ainsi les deux portions d'arc associés à la connecting line.
Une connecting line est supprimée si l'angle entre les orientations de deux portions d'arc est supérieur à _CL_EDGE_MAX_ANGLE_ ou bien si la distance de hausdorff entre elles est supérieure à _CL_EDGE_MAX_DIST_.
En outre, chaque connecting line isolée (non connectée à d'autres connecting lines) dont la longeur est inférieure à _CL_MIN_LENGTH_ et supprimée.

> Note : on désigne ici par orientation d'une polyligne, le vecteur défini par ses deux extrémitées.

![213](images/213_with_key.png)

#### 214 : UpdateGeomConnectingLines

##### Données de travail :

| table                          | entrée | sortie | entitée de travail | description                                                 |
|--------------------------------|--------|--------|--------------------|-------------------------------------------------------------|
| CL_TABLE                       |        | X      | X                  | table des connecting lines                                  |
| EDGE_TABLE_INIT                | X      |        |                    | table du réseau à traiter                                   |

##### Principaux opérateurs de calcul utilisés :
- app::calcul::CFeatGenerationOp

##### Description du traitement :

Paramètre utilisés: 
| paramètre                   | description                                                                                                               |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------|
| CL_SNAP_PROJ_CL_2_EDGE_DIST           | angle m
| LINKED_FEATURE_ID           | nom du champ indiquant les identifiants des arcs associés à une connecting line                                           |
EDGE_FICTITIOUS_NAME


calcul de la nouvelle géométrie des connecting lines par interpolation (géométrie moyenne) des deux portions d'arc associées. Sauf il un des deux arcs est fictifs et pas l'autre --> on prend la géom de l'arc fictif
La modification des géométries ayant créé des discontinuités on ré-établi les connexion légitimes  (cl connecté avant changement de géom et arcs associés identique ou connecté dans les pays respectifs)




TODO vérifier l'utilisation de EDGE_TABLE_INIT dans chacun des steps
TODO 
ign::geometry::LineString app::calcul::CFeatGenerationOp::_getGeomCL(
	epg::tools::MultiLineStringTool & mslBorder,
	ign::geometry::LineString const& lsToProject,
	double distMaxBorder,    -----------------------------------------------> pas utilisé
	double snapOnVertexBorder)

TODO ajouter CL_MIN_LENGTH