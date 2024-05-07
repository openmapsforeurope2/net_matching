#ifndef _APP_CALCUL_UTILS_ATTRIBUTEMERGER_H_
#define _APP_CALCUL_UTILS_ATTRIBUTEMERGER_H_


#include<ign/feature/Feature.h>

namespace app {
namespace calcul {
namespace utils {

	class AttributeMerger {

	public:

		AttributeMerger(std::string listAttr2concatName, std::string listAttrWName, std::string listAttrJsonName, std::string separtor);


		void addFeatAttributeMerger(
			ign::feature::Feature& featMerged,
			ign::feature::Feature& featAttrToAdd,
			std::string separator
		);


		std::set<std::string> getAttrNameToConcat() { return _sAttrNameToConcat;};
		std::set<std::string> getAttrNameW() { return _sAttrNameW;};
		std::set<std::string> getAttrNameJson() { return _sAttrNameJson;};


	private:

		std::set<std::string> _sAttrNameToConcat;
		std::set<std::string> _sAttrNameW;
		std::set<std::string> _sAttrNameJson;

	private:

		void _setListToSetAttr(std::string& listAttrName, std::set<std::string>& setAttrName,std::string separtor);


	};

}
}
}

#endif