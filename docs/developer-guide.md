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

Le processus de mise en cohérence des surfaces est décomposé en une succession d'étapes clés.
Afin d'orchestrer l'enchainement de ces étapes l'application utilise l'outil **epg::step::StepSuite** de la bibliothèque **libepg**. Ce dernier permet de lancer une succession de **epg::step::Step** dans lesquels sont implémentés les traitements de chaque étape.
Un code (numéro à trois chiffres) est attribué à chaque étape. Les étapes sont ordonnancées selon cette numérotation. Si une étape transforme les données sur lesquelles elle travaille, une ou plusieurs tables dédiées préfixées du code de l'étape sont créées. Ces créations sont réalisées en copiant les tables d'une étape antérieure (qui n'est pas nécessairement l'étape immédiatement antérieure, car toutes les étapes ne travaillent pas sur les mêmes données).
Ce fonctionnement permet de conserver les résultats intermédiaires du processus. Cela donne la possibilité d'arrêter et de reprendre le traitement en cours de processus et facilite le travail de d'analyse et de deboggage.


Les étapes qui composent le traitement de raccordement sont les suivantes :

201. initialisation du champ 'fictitious' des arcs du réseau
202. appariement des carrefours
204. généralisation des surfaces étroites en linéaires (les contours des surfaces sont constituées des arcs des 2 pays et les linéaires résultant de la fusion de ces arcs sont bi-nationaux)
210. génération des _'connectings lines'_ projetées aux frontières pour chacun des pays frontaliers
211. fusion des _'connecting lines'_ de chacun des pays en linéaires bi-nationaux
212. accrochage des _'connecting lines'_ dont les extrémités sont proches pour éviter les petites coupures
213. suppression des _'connecting lines'_ dont les couples de linéaires sont incohérents selon l'angle ou la distance
214. calcul de la géométrie des _'connecting lines'_ par interpolation 
220. connection du réseau aux _'connecting lines'_
230. import des _'connecting lines'_ dans le réseau
240. génération des _'connecting points'_
250. connection du réseau aux _'connecting points'_
255. généralisation des surfaces étroites en linéaires bi-nationaux
260. nettoyage des artefacts (faces étroites, antennes, arcs superposés...)
270. connection des antennes hors pays au réseau du pays voisin
280. nettoyage des artefacts (faces étroites, antennes, arcs superposés, arcs de petite taille...)

> _Précisions_ :
> - _'Connecting line' : arc résultant de la fusion de deux arcs de deux pays différents réprésentant la même portion de réseau._
> - _'Connecting point' : sommet représentant un point de passage du réseau à la frontière. Les réseaux des deux pays limitrophes doivent être connectés à ce point afin d'assurer la continuité topologique du référentiel européen._
> - _champ 'fictitious' : la valeur de ce champ est 'vraie' lorsqu'un arc est couvert par une géométrie surfacique représentant le même objet (certains réseaux possèdes deux représentations linéaire et surfacique modélisées par deux classes d'objets différentes. Tous les objets linéaires n'ont pas de représentation surfacique)._

L'outil **epg::step::StepSuite** donne la possibilité de ne lancer que certaines étapes ou une plage de plusieurs étapes.

La liste de l'ensemble des étapes qui constituent le processus de raccordement diffère selon la thématique traitée :
- hydrographie (hy) : étapes 201 à 280
- transport routier (tn) :  étapes 210 à 280
- transport ferré (ra) : étapes 240 à 280

## Configuration

L'outil s'appuie sur de nombreux paramètres de configuration permettant d'adapter le comportement des algorithmes en fonctions des spécificités nationales (sémantique, précision, échelle, conventions de modélisation...).

On trouve dans le [dossier de configuration](https://github.com/openmapsforeurope2/net_matching/tree/main/config) les fichiers suivants :

- epg_parameters.ini : regroupe des paramètres de base issus de la bibliothèque libepg qui constitue le socle de développement l'outil. Ce fichier est aussi le fichier chapeau qui pointe vers les autres fichiers de configurations.
- db_conf.ini : informations de connexion à la base de données.
- hy_theme_parameters.ini : configuration des paramètres spécifiques au thème hydrographie.
- ra_theme_parameters.ini : configuration des paramètres spécifiques au thème transport ferré.
- tn_theme_parameters.ini : configuration des paramètres spécifiques au thème transport routier.


## Lancement du traitement

L'outil s'utilise en ligne de commande.
Le traitement peut être lancé sur trois types de réseaux :
- routier (code tn)
- férré (code ra)
- hydrographique (code hy)

<br>

Les paramètres sont les suivants :
* c [obligatoire] : chemin vers le fichier de configuration
* T [obligatoire] : thème (doit être parmi les valeurs : tn, hy, ra)
* cc [obligatoire] : code pays double (exemple : be#fr)
* sp [obligatoire] : code de l'étape(s) à executer (exemples: 220 ou 220,240 ou 210-280...)

<br>

Exemple de lancement du traitement complet sur le couple de pays France (code pays 'fr') et Belgique (code pays 'be') pour le thème routier :
```
net_matching --c path/to/config/epg_parameters.ini --T tn --cc be#fr
```
A noter que l'on renseigne pour le paramètre --cc le code de la frontière séparant les deux pays à traiter. Le code pays est toujours composé de la même manière, c'est à dire en concaténant les codes <u>par ordre alphabétique</u>.

Exemple du lancement d'une seule étape :
```
net_matching --c path/to/config/epg_parameters.ini --cc be#fr --sp 240
```

Exemple de lancement d'une plage d'étapes :
```
net_matching --c path/to/config/epg_parameters.ini --cc be#fr --sp 240-260
```
Ici toutes les étapes de 240 à 260 (incluses) sont jouées.


## Les étapes - fonctionnement détaillé

### 201 : FillFictitiousField

Cette étape est spécifique au réseau hydrographique. Elle consiste à renseigner et/ou corriger le champ _EDGE_FICTITIOUS_NAME_ indiquant si un arc est fictif qui peut être manquant ou mal renseigné dans les donnnées sources livrées par les producteurs nationaux.
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

![201](images/201.png)


### 202 : JunctionMatching

Dans cette étape il est question d'appairer les carrefours (noeuds de degré supérieur à deux) des deux pays. Puis, pour chaque couple de carrefours appairés, réaliser les déplacements nécessaires pour leur localisation soit identique.
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

![202_a](images/202_a.png)

Afin d'éviter la création de rebroussement le principe constiste à projeter le point cible (point vers lequel l'extrémité de la polyligne doit être déplacée) sur la polyligne à déplacer, réaliser une découpe de la polyligne avec le point projeté, puis déplacer l'extrémité de la polyligne découpé vers le point cible.

![202_b](images/202_b.png)
<br>
![202_c](images/202_c.png)
