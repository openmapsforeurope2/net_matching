#ifndef _APP_GEOMETRY_TOOLS_LENGTHINDEXEDLINESTRING_H_
#define _APP_GEOMETRY_TOOLS_LENGTHINDEXEDLINESTRING_H_

#include <ign/ign_config.h>

#include <vector>

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
				/// \brief constructeur a partir d'une polyligne
				LengthIndexedLineString( ign::geometry::LineString const& lineString );
				/// \brief constructeur par recopie
				LengthIndexedLineString( LengthIndexedLineString const& other );
				/// \brief affectation
				LengthIndexedLineString& operator = ( LengthIndexedLineString const& other );
				/// \brief destructeur
				~LengthIndexedLineString();

				/// \brief renvoie la polyligne
				ign::geometry::LineString	const& 		getLineString() const ;
				/// \brief definit la polyligne a indexer
				void					setLineString( ign::geometry::LineString const& lineString );
				/// \brief renvoie l'abscisse d'un point
				double					getPointLocation( size_t numPoint ) const ;
				/// \brief renvoie une sous-partie de la linestring
				ign::geometry::LineString				getSubLineString( double const& sBegin, double const& sEnd );

				/// \brief renvoie les abscisses des points intermediaires
				std::vector< double > const & getPointAbscisses() const ;

				/// \brief renvoie la longueur de la polyligne
				double					length() const ;
				/// \brief renvoie le nombre de segments
				size_t					numLines() const ;

				/// \brief renvoie un point en fonction d'une abscisse curviligne
				ign::geometry::Point					locateAlong( double const& s ) const ;

				/// \brief renvoie l'abscisse curviligne du point le plus proche sur la
				/// polyligne du point p.
				double                  project( ign::geometry::Point const& p )const;

			private:
				/// \brief la linestring encapsulee
				ign::geometry::LineString				_lineString;
				/// \brief l'index
				std::vector< double >	_index;

                /// \brief gestion 3d si _lineString est 3D
                bool                    _is3d;
                
				/// \brief calcul l'index
				void	_computeIndex();

				size_t	_findLine( double const& s ) const ;

			};
		}
	}
}


#endif
