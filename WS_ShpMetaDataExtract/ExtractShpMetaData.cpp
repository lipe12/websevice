#include <iostream>
#include <strstream>
#include <string>
#include <vector>
#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "ogr_spatialref.h"
#include "tinyxml.h"
#include "tinystr.h"

using namespace std;

string ConvertToString(double val)
{  
	char dest[30];
	sprintf(dest, "%lf", val);
	string result(dest);
	return result;
}

string ConvertToString(int val)
{  
	char dest[30];
	sprintf(dest, "%d", val);
	string result(dest);
	return result;
}

int compare(const void *pFirst, const void *pSecond)
{
	int nA = *(( float *)pFirst);
	int nB = *(( float *)pSecond);
	return (nA<nB ? -1 : nA>nB ? 1 : 0); // 由小到大排序
	//return (nA<nB ? 1 : nA>nB ? -1 : 0); // 由大到小排序
}

TiXmlDocument* WriteMedtaDataToXML(vector<string> names, vector<string> values)
{
	//char *namespaceURI = "http://whu.edu.cn/ws/ExtractRasterMetaData";
	TiXmlDocument *doc = new TiXmlDocument();
	TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "UTF-8", "");
	TiXmlElement *root = new TiXmlElement("root");
	printf("%s",values[2].c_str());
	cout << values[2].c_str() <<endl;
	for (int i = 0; i < values.size(); i++ )
    {
		TiXmlElement *element = new TiXmlElement(names[i].c_str());
		TiXmlText *text = new TiXmlText(values[i].c_str());
		element->LinkEndChild(text);
		root->LinkEndChild(element);
	}
	doc->LinkEndChild(decl);
	doc->LinkEndChild(root);
	return doc;
}

TiXmlDocument* ExtractShpMetaData(char *shpFile)
{
	string ShpFileName;
	string ShpEPSGCODE;
	string ShpProjString;
	string ShpExtent;
	char *driverName = "ESRI Shapefile";
	OGRRegisterAll();
	OGRSFDriver *myDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
	OGRSpatialReference *srs = new OGRSpatialReference("");
	OGREnvelope *extent = new OGREnvelope();
	if (myDriver == NULL)
	{
		return NULL;
	}
	OGRDataSource *ds = myDriver->Open(shpFile, 0);
	if(ds == NULL)
	{
		return NULL;
	}
	OGRLayer *layer = ds->GetLayer(0);
	srs = layer->GetSpatialRef();
	int i = layer->GetExtent(extent, 1);
	OGRDataSource::DestroyDataSource(ds);

	// FILE NAME
    string fullname = (string)shpFile;
	int pos = fullname.find_last_of('/');
	string filename(fullname.substr(pos + 1));
	ShpFileName = filename;

	// EPSG CODE or PROJ STRING
	char *p4;
	//int result = srs->exportToProj4(&p4);
	int result = srs->exportToMICoordSys(&p4);
	string tmp = p4;
	cout << tmp << endl;
	const char *PROJCS = srs->GetAuthorityName("PROJCS");
	int wkid = -1;
	if (PROJCS != NULL && PROJCS == "EPSG")
	{
		sscanf(srs->GetAuthorityCode("PROJCS"),"%d", &wkid);
	}
	else if (srs->IsGeographic() == 1)
	{
		const char *GEOGCS = srs->GetAuthorityName("GEOGCS");
		if (GEOGCS != NULL && GEOGCS =="EPSG")
		{
			sscanf(srs->GetAuthorityCode("GEOGCS"),"%d", &wkid);
		}
	}
	strstream ss;
    string s;
    ss << wkid;
    ss >> s;
	ShpEPSGCODE = s;
    ShpProjString = p4;
	cout << "a" << endl;
	cout << p4 << endl;

	// EXTENT
	string left = ConvertToString(extent->MinX);
	string top = ConvertToString(extent->MaxY);
	string buttom = ConvertToString(extent->MinY);
	string right = ConvertToString(extent->MaxX);
	ShpExtent = left + " " + buttom + " " + right + " " + top;

	vector<string> ShpMetaData;
	ShpMetaData.push_back(ShpFileName);
	ShpMetaData.push_back(ShpEPSGCODE);
	ShpMetaData.push_back(ShpProjString);
	ShpMetaData.push_back(ShpExtent);

	vector<string> Names;
	Names.push_back("Filename");
	Names.push_back("EPSGCODE");
	Names.push_back("ProjString");
	Names.push_back("Exent");     

	return WriteMedtaDataToXML(Names, ShpMetaData);
}

int main(int argc, char *argv[])
{
	
	char *shpFile = argv[1];
	TiXmlDocument *doc2 = ExtractShpMetaData(shpFile);
	doc2->SaveFile(argv[2]);
	cout << "abo" <<endl;
	return 0;
}