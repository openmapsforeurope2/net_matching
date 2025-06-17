#ifndef _APP_CALCUL_CONNECTIVITYCORRECTOROP_H_
#define _APP_CALCUL_CONNECTIVITYCORRECTOROP_H_

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>

// SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>


namespace app{
namespace calcul{

	//--
	enum ENDING{
		START,
		END
	};

	//--
	struct PointCompare {
		bool operator()(ign::geometry::Point const& a, ign::geometry::Point const& b) const {
			return a.x() < b.x() || a.x() == b.x() && a.y() < b.y();
		}
	};

	/// @brief Opérateur pour améliorer la connectivité des réseaux nationaux
	class ConnectivityCorrectorOp {
		
	public:

		/// @brief Correction de la connectivité des réseaux
		/// @param verbose Code de la frontière traitée
		/// @param verbose Mode verbeux
		static void Compute(
			std::string const& borderCode,
			bool verbose
		);

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		bool                                               _verbose;

	private:

		//--
		ConnectivityCorrectorOp(bool verbose = false);

		//--
		~ConnectivityCorrectorOp();

		//--
		void _init();

		//--
		void _compute(
			std::string const& countryCode
		) const;

		//--
		void _createMissingConnections(
            std::string const& countryCode
        ) const;

		//--
		void _connectAtEndings(
            std::string const& countryCode
        ) const;

		//--
		ign::geometry::Point _getClosest(
            ign::geometry::Point const& pt, 
            std::vector<ign::geometry::Point> const& vPt
        ) const;

    };
}
}

#endif