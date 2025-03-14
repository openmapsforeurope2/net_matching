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

	/// \warning edges must have been displaced before to use this function
	/// @brief 
	/// @tparam GraphType 
	/// @param graph 
	/// @param v 
	/// @param vect 
	/// @param mOldNewEdges 
	/// @param sEdges2remove 
	/// @param sVertices2remove 
	/// @param withMerging 
	/// @param precision 
	/// @return 
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
