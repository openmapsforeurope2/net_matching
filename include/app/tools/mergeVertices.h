#ifndef _APP_TOOLS_MERGEVERTICES_H_
#define _APP_TOOLS_MERGEVERTICES_H_

// SOCLE
#include <ign/tools/stringtools.h>

// EPG
#include <epg/Context.h>
#include <epg/graph/concept/geom.h>
#include <epg/graph/concept/vertexFictitious.h>
#include <epg/graph/concept/countryCode.h>
#include <epg/graph/concept/setAttribute.h>
#include <epg/tools/StringTools.h>
#include <epg/graph/concept/feature.h>


namespace app{
namespace tools{

	/// \warning edges must have been displaced before to use this function
	/// @brief 
	/// @tparam GraphType 
	/// @param graph 
	/// @param v 
	/// @param vRef 
	/// @param mOldNewEdges 
	/// @param sEdges2remove 
	/// @param sVertices2remove 
	template< typename GraphType >
	void mergeVertices( 
		GraphType& graph,
		typename GraphType::vertex_descriptor v,
		typename GraphType::vertex_descriptor vRef,
		std::map<typename GraphType::edge_descriptor, typename GraphType::edge_descriptor> & mOldNewEdges,
		std::set<typename GraphType::edge_descriptor> & sEdges2remove,
		std::set<typename GraphType::vertex_descriptor> & sVertices2remove
	)
	{
		epg::Context* context = epg::ContextS::getInstance();

		typedef typename GraphType::vertex_descriptor           vertex_descriptor;
		typedef typename GraphType::edge_descriptor             edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor    oriented_edge_descriptor;

		std::vector< edge_descriptor > vIncidentEdges = graph.incidentEdges( v );

		typename std::vector< edge_descriptor >::const_iterator eit;
		for( eit = vIncidentEdges.begin() ; eit != vIncidentEdges.end() ; ++eit )
		{
			if (sEdges2remove.find(*eit) != sEdges2remove.end()) continue;
			//
			vertex_descriptor vSource = graph.source( *eit );
			vertex_descriptor vTarget = graph.target( *eit );
			if( vSource == v ) 
			{
				vSource = vRef;
			}
			if( vTarget == v )
			{
				vTarget = vRef;
			}
			
			//copie des proprietes
			typename GraphType::edge_properties properties = graph[*eit];

            std::vector< ign::geometry::Point > vIntermediatePoints = graph.intermediatePoints(*eit);

			oriented_edge_descriptor oe = graph.addEdge( vSource, vTarget, vIntermediatePoints, properties );

			std::vector<std::string> vOrigins = graph.origins(*eit);
			for (std::vector<std::string>::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit) {
				std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(*vit);
				if (foundInducedEdges.first) {
					for (typename std::vector<oriented_edge_descriptor>::iterator vit2 = foundInducedEdges.second.begin() ; vit2 != foundInducedEdges.second.end() ; ++vit2) {
						if (vit2->descriptor == *eit) {
							vit2->descriptor = oe.descriptor;
						}
					}
					graph.setInducedEdges(*vit, foundInducedEdges.second);
					graph.addOrigin( oe.descriptor , *vit );
				}
			}

			mOldNewEdges.insert(std::make_pair(*eit, oe.descriptor));
			sEdges2remove.insert(*eit);
		}

		sVertices2remove.insert(v);
	}
}
}

#endif
