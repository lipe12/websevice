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
	
	TiXmlDocument *doc = new TiXmlDocument();
	TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "UTF-8", "");
	TiXmlElement *root = new TiXmlElement("root");
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

TiXmlDocument* ExtractRasterMetaData(char *rasterFile, char *mode, int numQuantile)
{
	string rasterFileName;
	string rasterEPSGCODE;
 	string rasterProjString;
	string rasterExtent;
	string rasterMax;
	string rasterMin;
	string rasterQuantileBreaks = "";
	string rasterUniqueValues = "";

	GDALAllRegister();
	GDALDataset *ds;
	ds = (GDALDataset*)GDALOpen(rasterFile, GA_ReadOnly);
	if(ds == NULL)
	{
		return false;
	}
	const char *proj = ds->GetProjectionRef();
	double *geoTransform = new double[6];
	ds->GetGeoTransform(geoTransform);
	GDALRasterBand *band = ds->GetRasterBand(1);
	int tempNum = -1;
	int *pbSuccess = &tempNum;
	double noData = band->GetNoDataValue(pbSuccess);
	if(*pbSuccess == 0) { noData = -9999; }
	int xSize = ds->GetRasterXSize();
	int ySize = ds->GetRasterYSize();
	int count = xSize * ySize;
	GDALDataType datatype = band->GetRasterDataType();
	float *dataBuf = new float[count];
	band->RasterIO(GF_Read, 0, 0, xSize, ySize, dataBuf, xSize, ySize, datatype, 0, 0);
	
	// file name
	string fullname = (string)rasterFile;
	int pos = fullname.find_last_of('/');
	string filename(fullname.substr(pos + 1));
	rasterFileName = filename;
	
	// EPSG CODE or PROJ STRING
	OGRSpatialReference *srs = new OGRSpatialReference("");
	string wktStr(proj);
	char *wkt = new char[strlen(proj)+1];
	strcpy(wkt, proj);
	srs->importFromWkt(&wkt);
	char *p4;
	int result = srs->exportToProj4(&p4);
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
	rasterEPSGCODE = s;
    rasterProjString = p4;

	// EXTENT
	string left = ConvertToString(geoTransform[0]);
	string top = ConvertToString(geoTransform[3]);
	string buttom = ConvertToString(geoTransform[3] - geoTransform[1] * ySize);
	string right = ConvertToString(geoTransform[0] + geoTransform[1] * xSize);
	rasterExtent = left + " " + buttom + " " + right + " " + top;

	// MIN MAX
	double min = 0, max = 0, mean = 0, stddev = 0;
	band->ComputeStatistics(false, &min, &max, &mean, &stddev, NULL, NULL);
	rasterMax = ConvertToString(max);
	rasterMin = ConvertToString(min);

	// QUANTILES
	if (mode == "Q")
	{            
		double *QuantileBreaks = new double[numQuantile - 1];
		int dataCount = 0;
		for(int i = 0; i < count; i++)
		{
			if(dataBuf[i] != noData) { dataCount++; }
		}
		float *dataNewBuf = new float[dataCount];
		int index = 0;
		for(int i = 0; i < count; i++)
		{
			if(dataBuf[i] != noData)
			{
				dataNewBuf[index] = dataBuf[i]; 
				index++;
			}
		}
		qsort(dataNewBuf, dataCount, sizeof(float), compare);
        int numDataInQuantile = (int)((dataCount) / numQuantile + 0.5);
		for (int i = 0; i < numQuantile - 1; i++)
		{
			QuantileBreaks[i] = dataNewBuf[(i + 1) * numDataInQuantile];
			rasterQuantileBreaks = rasterQuantileBreaks + " " + ConvertToString(QuantileBreaks[i]);
		}
		int len = rasterQuantileBreaks.length();
		rasterQuantileBreaks = rasterQuantileBreaks.substr(1, len-1);
	}

	// Unique Values
	if (mode == "U")
	{
		vector<float> UniqueValues;
		int dataCount = 0;
		for(int i = 0; i < count; i++)
		{
			if(dataBuf[i] != noData) { dataCount++; }
		}
		float *dataNewBuf = new float[dataCount];
		int index = 0;
		for(int i = 0; i < count; i++)
		{
			if(dataBuf[i] != noData)
			{
				dataNewBuf[index] = dataBuf[i]; 
				index++;
			}
		}
		qsort(dataNewBuf, dataCount, sizeof(float), compare);

		UniqueValues.push_back(dataNewBuf[0]);
		rasterUniqueValues = ConvertToString((int)dataNewBuf[0]);

		for (int i = 1; i < dataCount - 1; i++)
		{
			if (dataNewBuf[i] != UniqueValues[UniqueValues.size() - 1])
			{
				UniqueValues.push_back(dataNewBuf[i]);
				rasterUniqueValues = rasterUniqueValues + " " + ConvertToString((int)dataNewBuf[i]);
			}

		}
	}
	GDALClose((GDALDatasetH)ds);	// close dataset

	vector<string> RasterMetaData;
	RasterMetaData.push_back(rasterFileName);
	RasterMetaData.push_back(rasterEPSGCODE);
	RasterMetaData.push_back(rasterProjString);
	RasterMetaData.push_back(rasterExtent);
	RasterMetaData.push_back(rasterMax);
	RasterMetaData.push_back(rasterMin);
	if (mode == "Q")
	{
		RasterMetaData.push_back(rasterQuantileBreaks);
	}
	if (mode == "U")
	{
		RasterMetaData.push_back(rasterUniqueValues);
	}

	vector<string> Names;
	Names.push_back("Filename");
	Names.push_back("EPSGCODE");
	Names.push_back("ProjString");
	Names.push_back("Exent");
	Names.push_back("Min");
	Names.push_back("Max");
	if (mode == "Q")
	{
		Names.push_back("QuantileBreaks");
	}
	if (mode == "U")
	{
		Names.push_back("UniqueValues");
	}
	return WriteMedtaDataToXML(Names, RasterMetaData);
}

int main(int argc, char *argv[])
{
	char *rasterFile = argv[1];
	char *outxmlFile = argv[2];
	char *Mode = argv[3];
	TiXmlDocument *doc1 = ExtractRasterMetaData(rasterFile, Mode, 15);
	doc1->SaveFile(outxmlFile);
	return 0;
}