#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "gdal_priv.h"
#include "ogrsf_frmts.h" //for ogr
#include "gdal_alg.h"	 //for GDALPolygonize
#include "cpl_conv.h"	 //for CPLMalloc()
#include "ogr_spatialref.h"
using namespace std;
void main(int argc, char *argv[])
{
	// ******************* 数据准备--载入测试数据 ******************* //

	// 参数示例
	//argc = 5;
	//argv[1] = "D:/RasterCut/input/dem1.tif";
	//argv[2] = "D:/1.txt"
	GDALAllRegister();
	char *str_inputRasterFileNames = argv[1];
	char *outFilepath = argv[2];
	GDALDataset *ds;
	ds = (GDALDataset*)GDALOpen(str_inputRasterFileNames, GA_ReadOnly);
	const char *GetProjectionRef = ds->GetProjectionRef();
	
	OGRSpatialReference *srs = new OGRSpatialReference("");
	string wktStr(GetProjectionRef);
	char *wkt = new char[strlen(GetProjectionRef) + 1];
	strcpy(wkt, GetProjectionRef);
	srs->importFromWkt(&wkt);
	char *p4;
	int result = srs->exportToProj4(&p4);
	const char *projs = srs->GetAuthorityName("PROJCS");
	const char *wkid = srs->GetAuthorityCode("PROJCS");
	string s(wkid);
	string epsgCode = wkid;
	int epsg = atoi(epsgCode.c_str());
	double left = 0;
	double right = 0;
	double top = 0;
	double bottom = 0;
	double geoTransform[6];
	ds->GetGeoTransform(geoTransform);
	left = geoTransform[0];
	top = geoTransform[3];
	right = geoTransform[0] + geoTransform[1] * ds->GetRasterXSize() + geoTransform[2] * ds->GetRasterYSize();
	bottom = geoTransform[3] + geoTransform[4] * ds->GetRasterXSize() + geoTransform[5] * ds->GetRasterYSize();
	/*char *wgs84 = "EPSG:4326";
	OGRSpatialReference oSRS;
	oSRS.SetWellKnownGeogCS(epsgCode.c_str());
	OGRSpatialReference oSRD;
	oSRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oSRD);
	coordTrans->Transform(1, &left, &top);
	int reprojected2 = coordTrans->Transform(1, &right, &bottom);*/
	OGRSpatialReference  oUTM, *poLatLong;
	OGRCoordinateTransformation *poTransform;
	oUTM.importFromEPSG(epsg);
	//poLatLong = oUTM.CloneGeogCS();
	poLatLong = new OGRSpatialReference();
	poLatLong->SetWellKnownGeogCS("WGS84");
	poTransform = OGRCreateCoordinateTransformation( &oUTM, poLatLong );
	poTransform->Transform(1, &left, &top);
	poTransform->Transform(1, &right, &bottom);
	ofstream outfile;
	outfile.open(outFilepath);
	outfile << left << " " << top << " " << right << " "<< bottom << endl;
	outfile.close();

}