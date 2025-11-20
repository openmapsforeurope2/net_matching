# net_matching

## Context

Open Maps For Europe 2 est un projet qui a pour objectif de développer un nouveau processus de production dont la finalité est la construction d'un référentiel cartographique pan-européen à grande échelle (1:10 000).

L'élaboration de la chaîne de production a nécessité le développement d'un ensemble de composants logiciels qui constituent le projet [OME2](https://github.com/openmapsforeurope2/OME2).

## Description

Le présent outil est dédié aux raccordement transfrontaliers des réseaux (réseaux routier, ferré, hydrographique...).

Lorsqu'elle est lancée l'application traite un couple de pays frontaliers. Pour raccorder l'ensemble du réseau d'un pays le programme doit être lancé successivement sur ses différentes frontières (en considérant l'ensemble de ses pays limitrophes).


## Fonctionnement

Le programme ne manipule pas directement les données de production. Les données à traiter, localisées autour de la frontière, sont extraites dans une table de travail. A l'issu du traitement les données dérivées sont injectées dans la table source en remplacement des données initiales.

Le processus de raccordement est décomposé en plusieurs étapes. Un numéro est attribué à chaque étape. Une table de travail préfixée de ce numéro est délivrée en sortie de chaque étape. Chaque étape prend en données d'entrées les tables de travail générées lors d'étapes antérieures.

Voici la liste de l'ensemble des étapes dont dispose l'outil :

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


## Utilisation

L'outil s'utilise en ligne de commande.
Le traitement peut être lancé sur trois types de réseaux:
- routier (code tn)
- férré (code ra)
- hydrographique (code hy)

Paramètres:
* c [obligatoire] : chemin vers le fichier de configuration
* s [obligatoire] : suffix de la table de travail
* as [optionnel] : suffix des tables de travail des surfaces (à utiliser uniquement lors du traitement de la classe d'objets watercourse_link)
* t [obligatoire] : nom de la classe d'objet (doit être parmi les valeurs : road_link, railway_link, watercourse_link)
* sp [obligatoire] : code de l'étape(s) à executer (exemples: 220 ou 220,240 ou 210-280...)
* arguments libres [obligatoire] : codes des deux pays frontaliers

<br>

Exemple d'appel pour lancer successivement l'ensemble des étapes du thème hydrographie sur la frontière franco-belge :
~~~
bin/net_matching --c path/to/config/epg_parameters.ini --s 20251118 --as 20251113 --t watercourse_link be fr
~~~

Exemple d'appel pour ne lancer qu'une seule étape :
~~~
bin/net_matching --c path/to/config/epg_parameters.ini --s 20251118 --as 20251113 --t watercourse_link --sp 250 be fr
~~~