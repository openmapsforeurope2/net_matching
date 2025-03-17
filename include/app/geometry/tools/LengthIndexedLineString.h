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
				/// \brief Constructeur par defaut
				LengthIndexedLineString();
				
				/// @brief Constructeur a partir d'une polyligne
				/// @param lineString 
				LengthIndexedLineString( ign::geometry::LineString const& lineString );

				/// @brief Constructeur par recopie
				/// @param other 
				LengthIndexedLineString( LengthIndexedLineString const& other );

				/// @brief Affectation
				/// @param other 
				/// @return 
				LengthIndexedLineString& operator = ( LengthIndexedLineString const& other );

				/// @brief Destructeur
				~LengthIndexedLineString();

				/// @brief Renvoie la polyligne
				/// @return 
				ign::geometry::LineString const& getLineString() const ;
				
				/// @brief Définit la polyligne a indexer
				/// @param lineString 
				void setLineString( ign::geometry::LineString const& lineString );

				/// @brief Renvoie l'abscisse d'un point
				/// @param numPoint index du point
				/// @return 
				double getPointLocation( size_t numPoint ) const ;

				/// @brief Renvoie une sous-partie de la linestring
				/// @param sBegin Abscisse curviligne de début
				/// @param sEnd Abscisse curviligne de fin
				/// @return 
				ign::geometry::LineString getSubLineString( double const& sBegin, double const& sEnd );

				/// @brief Renvoie les abscisses des points intermediaires
				/// @return Abscisses des points intermédiaires
				std::vector< double > const & getPointAbscisses() const ;

				/// @brief Renvoie la longueur de la polyligne
				/// @return Longueur
				double length() const ;

				/// @brief Renvoie le nombre de segments
				/// @return Nombre de segments
				size_t numLines() const ;

				/// @brief Renvoie un point en fonction d'une abscisse curviligne
				/// @param s Abscisse curviligne
				/// @return Point
				ign::geometry::Point locateAlong( double const& s ) const ;

				/// @brief Renvoie l'abscisse curviligne du point le plus proche sur la polyligne du point p.
				/// @param p Point
				/// @return Abscisse
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
