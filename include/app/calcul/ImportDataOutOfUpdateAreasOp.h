#ifndef _APP_CALCUL_IMPORTDATAOUTOFUPDATEAREASOP_H_
#define _APP_CALCUL_IMPORTDATAOUTOFUPDATEAREASOP_H_

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>

// IGN
#include <ign/feature/sql/FeatureStorePostgis.h>


namespace app{
namespace calcul{

	/// @brief Classe utilisée pour l'appairage des carrefours
	class ImportDataOutOfUpdateAreasOp {

	public:

		/// @brief Constructeur
		/// @param verbose Mode verbeux
		ImportDataOutOfUpdateAreasOp(bool verbose = false);

		/// @brief Destructeur
		~ImportDataOutOfUpdateAreasOp();

		/// @brief Importe les données hors zone de mise à jour reservées lors de l'extraction initiale
		/// @param verbose Mode verbeux
		static void Compute(bool verbose = false);


	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsOfua;
		//--
		epg::log::EpgLogger*                               _logger;
		//--
		bool                                               _verbose;

		
	private:

		//--
		void _init();

		//--
		void _compute() const;

	};

}
}

#endif