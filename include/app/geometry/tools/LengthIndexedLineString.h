#ifndef _APP_GEOMETRY_TOOLS_LENGTHINDEXEDLINESTRING_H_
#define _APP_GEOMETRY_TOOLS_LENGTHINDEXEDLINESTRING_H_

#include <vector>

// SOCLE
#include <ign/ign_config.h>
#include <ign/geometry/LineString.h>


namespace app {
	namespace geometry {
		namespace tools {
			/// \brief ign::geometry::LineString indexee sur la longueur permettant d'interpoler les points
			/// en fonction de leurs abscisse
			class LengthIndexedLineString {

			public:
				/// \brief constructeur par defaut
				LengthIndexedLineString();
				
				/// @brief constructeur a partir d'une polyligne
				/// @param lineString 
				LengthIndexedLineString( ign::geometry::LineString const& lineString );

				/// @brief constructeur par recopie
				/// @param other 
				LengthIndexedLineString( LengthIndexedLineString const& other );

				/// @brief affectation
				/// @param other 
				/// @return 
				LengthIndexedLineString& operator = ( LengthIndexedLineString const& other );

				/// @brief destructeur
				~LengthIndexedLineString();

				/// @brief renvoie la polyligne
				/// @return 
				ign::geometry::LineString const& getLineString() const ;
				
				/// @brief definit la polyligne a indexer
				/// @param lineString 
				void setLineString( ign::geometry::LineString const& lineString );

				/// @brief renvoie l'abscisse d'un point
				/// @param numPoint 
				/// @return 
				double getPointLocation( size_t numPoint ) const ;

				/// @brief renvoie une sous-partie de la linestring
				/// @param sBegin 
				/// @param sEnd 
				/// @return 
				ign::geometry::LineString getSubLineString( double const& sBegin, double const& sEnd );

				/// @brief renvoie les abscisses des points intermediaires
				/// @return 
				std::vector< double > const & getPointAbscisses() const ;

				/// @brief renvoie la longueur de la polyligne
				/// @return 
				double length() const ;

				/// @brief renvoie le nombre de segments
				/// @return 
				size_t numLines() const ;

				/// @brief renvoie un point en fonction d'une abscisse curviligne
				/// @param s 
				/// @return 
				ign::geometry::Point locateAlong( double const& s ) const ;

				/// @brief renvoie l'abscisse curviligne du point le plus proche sur la polyligne du point p.
				/// @param p 
				/// @return 
				double project( ign::geometry::Point const& p )const;

			private:
				//--
				ign::geometry::LineString  _lineString;
				//--
				std::vector< double >      _index;
                //--
                bool                       _is3d;

			private:
				//--
				void _computeIndex();
				
				//--
				size_t _findLine( double const& s ) const ;

			};
		}
	}
}


#endif
