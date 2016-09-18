#include "ogrsf_frmts.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h" 
#include <string>
#include <iostream>
#include <strstream>
#include <exception>
using namespace std;

void convertToShp(double longitude, double latitude, char *outshp)
{
	
   // double *TransferedLongLat = transfer(longitude, latitude);	
	
	//使属性表字段支持中文
	CPLSetConfigOption("SHAPE_ENCODING","");
	OGRRegisterAll();//注册所有的驱动
	//创建ESRI shp文件
	char *pszDriverName = "ESRI Shapefile";
	//调用对Shape文件读写的Driver
	OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName);
	if (poDriver == NULL)
	{
		cout<<pszDriverName<<"驱动不可用！"<<endl;
		return;
	}
	//创建数据源
	OGRDataSource *poDs = poDriver->CreateDataSource(outshp, NULL);
	if (poDs == NULL)
	{
		cout<<"DataSource Creation Error"<<endl;
		return;
	}
	//创建图层Layer
	string outShapName = outshp;
	string layerName = outShapName.substr(0, outShapName.length()-4);
	//layerName.c_str()表示将string转为char *类型
	//参数说明：新图层名称，坐标系，图层的几何类型，创建选项，与驱动有关
	OGRLayer *poLayer = poDs->CreateLayer(layerName.c_str(), NULL, wkbPoint, NULL);
	if (poLayer == NULL)
	{
		cout<<"Layer Creation Failed"<<endl;
		OGRDataSource::DestroyDataSource(poDs);
		return;
	}
	//下面创建属性表，我们在属性表中创建两列数据即可
	//先创建一个“ID”整型属性
	OGRFieldDefn oFieldId("ID", OFTInteger);
	oFieldId.SetWidth(10);
	poLayer->CreateField(&oFieldId);
	//name
	OGRFieldDefn oFieldname("Dist_moved", OFTString);
	oFieldId.SetWidth(50);
	poLayer->CreateField(&oFieldname);
	//再创建一个"X"double属性
	OGRFieldDefn oFieldX("X", OFTString);
	oFieldX.SetWidth(50);
	poLayer->CreateField(&oFieldX);
	//再创建一个"Y"double属性
	OGRFieldDefn oFieldY("Y", OFTString);//OFTReal
	oFieldY.SetWidth(50);
	poLayer->CreateField(&oFieldY);
	//创建一个feature
	OGRFeature *poFeature; 	
	poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());//GetLayerDefn()获取当前图层的属性表结构
	//给属性表中我们刚创建的列赋值
	int i = 0;
	poFeature->SetField("ID", i);
	poFeature->SetField("Dist_moved", i);
	poFeature->SetField("X", longitude);
	poFeature->SetField("Y", latitude);
	i++;
	//创建一个OGRPoint对象
	OGRPoint point;
	point.setX(longitude);
	point.setY(latitude);
	//point.setZ(0);
	
	poFeature->SetGeometry(&point);

	if(poLayer->CreateFeature(poFeature) != OGRERR_NONE )
	{
		printf( "Failed to create feature in shapefile.\n" );
		exit( 1 );
	}
	OGRFeature::DestroyFeature(poFeature);
	OGRDataSource::DestroyDataSource(poDs);
	
}
//======= 经纬度转化为投影坐标=============
double* transfer(double longitude, double latitude)
{
	OGRSpatialReference oSourceSRS;
	//EPSG code 代表特定的椭球体、单位、地理坐标系或投影坐标系等信息
	//This method will initialize the spatial reference based on the passed in EPSG GCS or PCS code
	oSourceSRS.importFromEPSG(4326);//EPSG:4326代表地理坐标系WGS1984
	OGRSpatialReference oTargetSRS;
	oTargetSRS.importFromEPSG(2029);
	OGRCoordinateTransformation *poTransform;
	poTransform = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
	if (poTransform == NULL)
	{
		cout<<"poTransform is null"<<endl;
		exit(1);
	}	
	if (!poTransform->Transform(1, &longitude, &latitude, NULL))
	{
		cout<<"transform failed"<<endl;
		exit(1);
	}
	//poTransform->Transform(1, &longitude, &latitude, NULL);
	double *inout = new double[2];
	inout[0] = longitude;
	inout[1] = latitude;
	return inout;
}
//这个函数和上面的transfer函数的功能是一样的
double* transfer2(double longitude, double latitude)
{
	OGRSpatialReference oSourceSRS;	
	oSourceSRS.SetWellKnownGeogCS( "WGS84" );
	OGRSpatialReference oTargetSRS;	
	oTargetSRS.SetWellKnownGeogCS("WGS84");
	oTargetSRS.SetUTM(17);	
	OGRCoordinateTransformation *poTransform;
	poTransform = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
	if (poTransform == NULL)
	{
		cout<<"poTransform is null"<<endl;
		exit(1);
	}	
	if (!poTransform->Transform(1, &longitude, &latitude))
	{
		cout<<"transform failed"<<endl;
		exit(1);
	}
	//poTransform->Transform(1, &longitude, &latitude, NULL);
	double *inout = new double[2];
	inout[0] = longitude;
	inout[1] = latitude;
	return inout;
}
int main(int argc, char *argv[])
{		
	char *xCoordinate = argv[1];
	char *yCoordinate = argv[2];
	char *outShp = argv[3];
	double x = atof(xCoordinate);
	double y = atof(yCoordinate);
	convertToShp(x, y, outShp);
	cout<<"success! file is saved to "<<outShp<<endl;

	/*
	double *transferedLongLat = transfer(116.246742, 40.022211);
	double *transferedLongLat2 = transfer2(116.246742, 40.022211);
	cout<<"转换后的投影坐标为："<<transferedLongLat[0]<<","<<transferedLongLat[1]<<endl;	
	cout<<"使用第二个函数转换："<<transferedLongLat2[0]<<","<<transferedLongLat2[1]<<endl;	
	*/
	/*double x = 735847.112853;
	double y = 3428713.99719;
	const char *outShp = "E:\\Exercise\\test\\convertResult\\outshp.shp";*/			
	
	getchar();
}