# net_matching

## Context: Ome2

Thèmes

## Description

L'outil net_matching est dédié au raccordements transfrontaliers des réseaux

## Configuration


## Etapes

Travail sur une frontière (2 pays)
Etape préliminaire d'extraction du réseau autour des frontières

On travail sur la table de travail. A chaque étape une table de travail intermédiaire préfixée du numéro d'étape est créée.

201 : initialisation du champ 'fictitious' des arcs du réseau
202 : apparairage des carrefour
204 : généralisation des surfaces étroites en linéaire (les surfaces sont constitués des arcs de 2 pays, les linéaires résultant sont des 'connecting lines' bi-nationales)
210 : génération des 'connectings lines' projetées aux frontières
211 : fusion des 'connecting lines'
212 :
213 :
214 : calcul de la géométrie des 'connecting lines' par interpolation
220 : connection du réseau aux 'connecting lines'
230 : import des 'connecting lines' dans le réseau
240 : génération des 'connecting points'
250 : connection du réseau aux 'connecting points'
255 : généralisation des surfaces étroites en linéaire (les surfaces sont constitués des arcs de 2 pays, les linéaires résultant sont des 'connecting lines' bi-nationales)
260 : nettoyage des artefacts (faces étroites, antennes, arcs superposés...)
270 : connection des antennes hors pays au réseau du pays voisin
280 : nettoyage des artefacts (faces étroites, antennes, arcs superposés, arcs de petite taille...)


Voici quelles sont les étapes jouées en fonction du thème traité:
- hydrographie (hy) : étapes 201 à 280
- transport routier (tn) :  étapes 210 à 280
- transport ferré (ra) : étapes 240 à 280


## dépendances
## compilation
## utilisation

Paramètres:
* c [obligatoire] : chemin vers le fichier de configuration
* T [obligatoire] : thème (doit être parmi les valeurs : tn, hy, ra)
* cc [obligatoire] : code pays double (exemple : be#fr)
* sp [obligatoire] : étape(s) à executer (exemples: 220 ; 220,240 ; 210-280)

Exemples d'appels:

~~~
bin/net_matching --c path/to/config/epg_parmaters.ini --cc be#fr --T hy
bin/net_matching --c path/to/config/epg_parmaters.ini --cc be#fr --T hy --sp 250
~~~

## Définition

'Connecting line' : 
'Connecting point' :
champ 'fictitious' :

