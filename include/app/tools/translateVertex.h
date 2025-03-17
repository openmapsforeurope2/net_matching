#ifndef _APP_TOOLS_TRANSLATEVERTEX_H_
#define _APP_TOOLS_TRANSLATEVERTEX_H_

// APP
#include <app/tools/mergeVertices.h>

// EPG
#include <epg/Context.h>
#include <epg/graph/concept/geom.h>
#include <epg/graph/concept/vertexFictitious.h>
#include <epg/graph/concept/countryCode.h>
#include <epg/tools/StringTools.h>

namespace app{
namespace tools{

	/// @brief Déplace le sommet d'un graph et réalise la fusion si sa nouvelle position
	/// correspond à la position (à la précision près) d'un sommet existant.
	/// Note : les arcs doivent avoir été préalablement déplacés avant d'utiliser cette fonction.
	/// @tparam GraphType 
	/// @param graph Graph à modifier
	/// @param v Identifiant du sommet à déplacer
	/// @param vect Vecteur de déplacement
	/// @param mOldNewEdges Mapping entre les anciens ars (à supprimer) et les nouveaux.
	/// @param sEdges2remove Liste des arcs à supprimer
	/// @param sVertices2remove Liste des sommets à supprimer
	/// @param withMerging Booléen indiquant si l'éventuelle fusion doit être réalisée
	/// @param precision Précision
	/// @return Pair indiquant si une fusion a été réalisée, et, si oui, l'identifiant du sommet 
	/// auquel le sommet déplacé a été fusionné.
	template< typename GraphType >
	std::pair< bool/*merged*/, typename GraphType::vertex_descriptor > translateVertex( 
		GraphType& graph,
		typename GraphType::vertex_descriptor v,
		ign::math::Vec2d const& vect,
		std::map<typename GraphType::edge_descriptor, typename GraphType::edge_descriptor> & mOldNewEdges,
		std::set<typename GraphType::edge_descriptor> & sEdges2remove,
		std::set<typename GraphType::vertex_descriptor> & sVertices2remove,
		bool withMerging,
		double precision = 1e-7
	)
	{
		typedef typename GraphType::vertex_descriptor  vertex_descriptor;
		typedef typename GraphType::edge_descriptor    edge_descriptor;

		epg::Context* context = epg::ContextS::getInstance();
		ign::geometry::Point oldGeom = graph.getGeometry(v);
		ign::geometry::Point newGeom( oldGeom.x()+vect.x(), oldGeom.y()+vect.y(), oldGeom.z() );

		//on recherche une eventuelle fusion
		std::set< vertex_descriptor > sVertices;
		graph.getVertices( newGeom.getEnvelope().expandBy( precision ), sVertices );
		
		double distMin = precision;
		bool merged = false;
		vertex_descriptor minVertex = GraphType::nullVertex();
		typename std::set< vertex_descriptor >::const_iterator vit;
		for( vit = sVertices.begin() ; vit != sVertices.end() ; ++vit )
		{
			if( *vit == v ) continue;
			if( sVertices2remove.find(*vit) != sVertices2remove.end()) continue;
			double distance = newGeom.distance(graph.getGeometry(*vit));
			if( distance < distMin )
			{
				distMin = distance;
				minVertex = *vit;
				merged = true;
			}
		}
		if( merged && withMerging )
		{
			mergeVertices( graph, v, minVertex, mOldNewEdges, sEdges2remove, sVertices2remove );
		}
		else
		{
			graph.setGeometry(v, newGeom );
		}

		return std::make_pair( merged, minVertex );
	}

}
}

#endif
