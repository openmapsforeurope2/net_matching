#include <app/calcul/utils/attributeMerger.h>


#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>

#include <epg/tools/StringTools.h>

//DEBUG
#include <epg/log/EpgLogger.h>


app::calcul::utils::AttributeMerger::AttributeMerger(std::string listAttr2concatName, std::string listAttrWName, std::string listAttrJsonName, std::string separtor )
{

	_setListToSetAttr(listAttr2concatName, _sAttrNameToConcat, separtor);

	//DEBUG
	epg::log::EpgLogger* logger = epg::log::EpgLoggerS::getInstance();
	logger->log(epg::log::DEBUG, "C1");
	for (std::set<std::string>::const_iterator sit = _sAttrNameW.begin() ; sit != _sAttrNameW.end() ; ++sit ) {
		logger->log(epg::log::DEBUG, *sit);
	}

	_setListToSetAttr(listAttrWName, _sAttrNameW, separtor);

	//DEBUG
	logger->log(epg::log::DEBUG, "C2");
	for (std::set<std::string>::const_iterator sit = _sAttrNameW.begin() ; sit != _sAttrNameW.end() ; ++sit ) {
		logger->log(epg::log::DEBUG, *sit);
	}

	_setListToSetAttr(listAttrJsonName, _sAttrNameJson, separtor);

}

app::calcul::utils::AttributeMerger::~AttributeMerger()
{
}


void app::calcul::utils::AttributeMerger::_setListToSetAttr(std::string& listAttrName, std::set<std::string>& setAttrName, std::string separtor)
{
	//DEBUG
	epg::log::EpgLogger* logger = epg::log::EpgLoggerS::getInstance();
	logger->log(epg::log::DEBUG, "B1");
	logger->log(epg::log::DEBUG, listAttrName);
	logger->log(epg::log::DEBUG, separtor);

	std::vector<std::string> vAttrName;
	epg::tools::StringTools::Split(listAttrName, separtor, vAttrName);
	for (size_t i = 0; i < vAttrName.size(); ++i){
		logger->log(epg::log::DEBUG, vAttrName[i]);

		setAttrName.insert(vAttrName[i]);
		logger->log(epg::log::DEBUG, "B2");
	}
		
	
}


void app::calcul::utils::AttributeMerger::addFeatAttributeMerger(
	ign::feature::Feature& featMerged,
	ign::feature::Feature& featAttrToAdd,
	std::string separator
)
{
	//DEBUG
	epg::log::EpgLogger* logger = epg::log::EpgLoggerS::getInstance();
	logger->log(epg::log::DEBUG, "A1");
	logger->log(epg::log::DEBUG, featMerged.getId());
	logger->log(epg::log::DEBUG, featAttrToAdd.getId());

	ign::feature::FeatureType featTypMerged = featMerged.getFeatureType();
	//ign::feature::FeatureType featTyp2 = feat2.getFeatureType();
	//test si featTyp1 == featTyp2 sinon msg d'erreur
	std::vector<std::string> vAttrNames;
	featTypMerged.getAttributeNames(vAttrNames);

	//DEBUG
	logger->log(epg::log::DEBUG, "A2");

	for (size_t i = 0; i < vAttrNames.size(); ++i) {
		std::string attrValueMerged;
		std::string attrName = vAttrNames[i];

		//DEBUG
		logger->log(epg::log::DEBUG, attrName);
		for (std::set<std::string>::const_iterator sit = _sAttrNameW.begin() ; sit != _sAttrNameW.end() ; ++sit ) {
			logger->log(epg::log::DEBUG, *sit);
		}

		if (_sAttrNameW.find(attrName) != _sAttrNameW.end()) //on ne fusionne pas les attributs de travail
			continue;

		//DEBUG
		logger->log(epg::log::DEBUG, "A3");

		if (_sAttrNameJson.find(attrName) != _sAttrNameJson.end()) { // on ajoute les json aux jsonArray

			//DEBUG
			logger->log(epg::log::DEBUG, "A4");

			std::string attrValueToMerge = featMerged.getAttribute(attrName).toString();
			std::string attrValueToAdd = featAttrToAdd.getAttribute(attrName).toString();
			QJsonDocument attrValueToMergeJsonDoc = QJsonDocument::fromJson(QString(attrValueToMerge.c_str()).toUtf8());
			QJsonDocument attrValueToAddJsonDoc = QJsonDocument::fromJson(QString(attrValueToAdd.c_str()).toUtf8());

			//DEBUG
			logger->log(epg::log::DEBUG, "A5");

			if (attrValueToMerge == "{}")
				attrValueMerged = attrValueToAdd;
			else if (attrValueToAdd == "{}")
				attrValueMerged = attrValueToMerge;
			else {

				//DEBUG
				logger->log(epg::log::DEBUG, "A6");

				QJsonArray attrArrayJsonToMerge;
				QJsonArray attrArrayJsonToAdd;

				if (attrValueToMergeJsonDoc.isArray())
					attrArrayJsonToMerge = attrValueToMergeJsonDoc.array();
				if (attrValueToAddJsonDoc.isArray())
					attrArrayJsonToAdd = attrValueToAddJsonDoc.array();

				QJsonArray attrValueMergedJsonArray;

				for (QJsonArray::iterator itA = attrArrayJsonToMerge.begin(); itA != attrArrayJsonToMerge.end(); ++itA)
					attrValueMergedJsonArray.push_back(*itA);
				for (QJsonArray::iterator itA = attrArrayJsonToAdd.begin(); itA != attrArrayJsonToAdd.end(); ++itA)
					attrValueMergedJsonArray.push_back(*itA);

				//DEBUG
				logger->log(epg::log::DEBUG, "A6");

				QJsonDocument attrValueMergedJsonDoc;
				attrValueMergedJsonDoc.setArray(attrValueMergedJsonArray);
				attrValueMerged = attrValueMergedJsonDoc.toJson(QJsonDocument::Compact).toStdString().c_str();

				//DEBUG
				logger->log(epg::log::DEBUG, "A7");
			}
			//DEBUG
			logger->log(epg::log::DEBUG, "A8");
			featMerged.setAttribute(attrName, ign::data::String(attrValueMerged));
			//DEBUG
			logger->log(epg::log::DEBUG, "A9");
		}
		else {
			//DEBUG
			logger->log(epg::log::DEBUG, "A11");

			std::string attrValueToMerge = featMerged.getAttribute(attrName).toString();
			std::string attrValueToAdd = featAttrToAdd.getAttribute(attrName).toString();
			if (attrValueToMerge == "void_unk" || attrValueToMerge == "-32768")
				attrValueMerged = attrValueToAdd;
			else if (attrValueToAdd == "void_unk" || attrValueToAdd == "-32768")
				attrValueMerged = attrValueToMerge;
			else if (attrValueToMerge != attrValueToAdd)
				attrValueMerged = attrValueToMerge + separator + attrValueToAdd;
			else
				attrValueMerged = attrValueToMerge;

			//DEBUG
			logger->log(epg::log::DEBUG, "A12");
			featMerged.setAttribute(attrName, ign::data::String(attrValueMerged));
			//DEBUG
			logger->log(epg::log::DEBUG, "A13");
		}

	}
}
