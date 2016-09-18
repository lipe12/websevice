#include "RasterLayer.h"
#include "gdal_priv.h"
#include "ogr_spatialref.h"
#include <iostream>
#include "Utility.h"

using namespace std;

CRasterLayer::CRasterLayer()
{
	//m_iLayerID = -1;
	m_strLayerName = "";
	m_strFileName = "";
	m_pDataBuf = NULL;
	NumberOfRows = 0;
	NumberOfCols = 0;
	m_fNullDataValue = NOFLOATVALUE;

}

//copy construction function
CRasterLayer::CRasterLayer(CRasterLayer &anEnvLayer) 
{
	m_pDataBuf = NULL;
	//m_iLayerID = anEnvLayer.getLayerID();
	m_strLayerName = anEnvLayer.getLayerName();
	m_strFileName = anEnvLayer.getFileName();

	if(anEnvLayer.m_pDataBuf != NULL)
		this->readData();
}



CRasterLayer::CRasterLayer(string strLayerName, string strFileName)
{
	m_pDataBuf = NULL;
	m_strLayerName = strLayerName;
	m_strFileName = strFileName;
	this->readHeader();
}


CRasterLayer::CRasterLayer(string strFileName)
{
	m_pDataBuf = NULL;
	m_strLayerName = extract_filename((char *)m_strFileName.data());
	m_strFileName = strFileName;
	this->readHeader();
}



CRasterLayer::~CRasterLayer()
{
	if(NULL != m_pDataBuf)
	{
		delete []m_pDataBuf;
		m_pDataBuf = NULL;
	}
}


void CRasterLayer::setLayerName(string strLayerName)
{
	m_strLayerName = strLayerName;
}

string CRasterLayer::getLayerName()
{
	return m_strLayerName;
}

void CRasterLayer::setFileName(string strFileName)
{
	m_strFileName = strFileName;
}

string CRasterLayer::getFileName()
{
	return m_strFileName;
}

int CRasterLayer::getNumberOfRows()
{
	return NumberOfRows;

}

int CRasterLayer::getNumberOfColumns()
{

	return NumberOfCols;

}

float CRasterLayer::getCellSize()
{
	return CellSize;
}


float CRasterLayer::getData(int iRow, int iCol)
{
	if(m_pDataBuf != NULL)
	{
		if(fabs(m_pDataBuf[iRow * NumberOfCols +iCol] - m_fNullDataValue) > VERYSMALL)
			return m_pDataBuf[iRow * NumberOfCols +iCol];
	}
	return NOFLOATVALUE;

}

float CRasterLayer::getData(double x,double y)
{
	if(m_pDataBuf != NULL)
	{
		int iCol = (x - Xmin)/CellSize;
		int iRow = (Ymax - y)/CellSize;

		if (iCol >= 0 && iCol < NumberOfCols && iRow >= 0 && iRow < NumberOfRows)
		{
			if(fabs(m_pDataBuf[iRow * NumberOfCols +iCol] - m_fNullDataValue) > VERYSMALL)
				return m_pDataBuf[iRow * NumberOfCols +iCol];
		}
	}
	return NOFLOATVALUE;

}


int CRasterLayer::readHeader()
{

	GDALAllRegister();


	GDALDataset *poDataset;

	poDataset = (GDALDataset *) GDALOpen(m_strFileName.data(), GA_ReadOnly);
	if(poDataset != NULL)
	{

		long bandNum;
		string DataType;
		double AdfGeoTransform[6];

		int CurrentBand;

		NumberOfCols = poDataset->GetRasterXSize(); 
		NumberOfRows = poDataset->GetRasterYSize();
		bandNum=poDataset->GetRasterCount();

		if(bandNum > 1)
		{
			cout <<"There are more than one band. Only the first one will be converted." <<endl;
		}


		if( poDataset->GetProjectionRef()  != NULL )
		{
			SpatialRef =poDataset->GetProjectionRef(); 
		}

		if( poDataset->GetGeoTransform(AdfGeoTransform ) == CE_None )
		{		
			Xmin = AdfGeoTransform[0];
			Ymax = AdfGeoTransform[3] ;
			CellSize = AdfGeoTransform[1];
			Xmax = AdfGeoTransform[0] + NumberOfCols*AdfGeoTransform[1] + NumberOfRows*AdfGeoTransform[2];
			Ymin = AdfGeoTransform[3] + NumberOfCols*AdfGeoTransform[4] + NumberOfRows*AdfGeoTransform[5];

		}

		GDALRasterBand  *poBand=NULL;
		poBand = poDataset->GetRasterBand(1);	
		if(poBand)
		{
			DataType=GDALGetDataTypeName(poBand->GetRasterDataType());
			this->m_fNullDataValue = poBand->GetNoDataValue();
			string str = poBand->GetDescription();
			str = poBand->GetUnitType();
	
		}


		GDALClose( (GDALDatasetH) poDataset);
		return 0;
	}

	GDALClose( (GDALDatasetH) poDataset);
	return 1;
}


int CRasterLayer::releaseDataBuf()
{
	if(NULL != m_pDataBuf)
	{
		delete []m_pDataBuf;
		m_pDataBuf = NULL;
	}
	return 1;
}

int CRasterLayer::readData()
{
	if(NULL != m_pDataBuf)
	{
		cout << "The data have already been read in." << endl;
		return 1;
	}

	GDALDataset *poDataset;

	poDataset = (GDALDataset *) GDALOpen(m_strFileName.data(), GA_ReadOnly);
	if(poDataset != NULL)
	{

		long bandNum;

		string Projection ;
		string DataType;
		double AdfGeoTransform[6];

		int CurrentBand;

		NumberOfCols = poDataset->GetRasterXSize(); 
		NumberOfRows = poDataset->GetRasterYSize();
		bandNum=poDataset->GetRasterCount();

		if(bandNum > 1)
		{
			cout <<"There are more than one band. Only the first one will be converted." <<endl;
		}


		if( poDataset->GetProjectionRef()  != NULL )
		{
			Projection=poDataset->GetProjectionRef(); 
		}

		if( poDataset->GetGeoTransform( AdfGeoTransform ) == CE_None )
		{		
			Xmin = AdfGeoTransform[0];
			Ymax = AdfGeoTransform[3] ;
			CellSize = AdfGeoTransform[1];
			Xmax = AdfGeoTransform[0] + NumberOfCols*AdfGeoTransform[1] + NumberOfRows*AdfGeoTransform[2];
			Ymin = AdfGeoTransform[3] + NumberOfCols*AdfGeoTransform[4] + NumberOfRows*AdfGeoTransform[5];

		}

		GDALRasterBand  *poBand=NULL;


		//CurrentBand=1;
		poBand = poDataset->GetRasterBand(1);	
		if(poBand)
		{
			DataType=GDALGetDataTypeName(poBand->GetRasterDataType());
			this->m_fNullDataValue = poBand->GetNoDataValue();
			string str = poBand->GetDescription();
			str = poBand->GetUnitType();
			//int dataType=0;
			/*if(DataType=="Byte" || DataType == "B")
				this->m_iDataType = 0;
			else if(DataType == "Float32" || DataType == "f" || DataType == "Float64")
				this->m_iDataType = 1;
			else 
				this->m_iDataType = 2;*/

			m_pDataBuf = new float[NumberOfRows * NumberOfCols];
			if (CE_None==poBand->RasterIO(GF_Read, 0,0, NumberOfCols, NumberOfRows, m_pDataBuf, NumberOfCols, NumberOfRows, GDT_Float32, 0, 0 ))

			{
				GDALClose( (GDALDatasetH) poDataset);
				return 0;
			}
			else 
			{
				GDALClose( (GDALDatasetH) poDataset);
				return 1;
			}
		}
	}


	GDALClose( (GDALDatasetH) poDataset);
	return 1;



}

int CRasterLayer::wirteData(string format)
{
	const char *pszFormat = format.data();
	GDALDriver *poDriver;
	char **papszMetadata;

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	if( poDriver == NULL )
		exit( 1 );



	GDALDataset *poDstDS;       
	char **papszOptions = NULL;

	poDstDS = poDriver->Create( this->m_strFileName.data(), NumberOfCols, NumberOfRows, 1, GDT_Float32, 
		papszOptions );


	double adfGeoTransform[6] = {Xmin, CellSize, 0, Ymax, 0, -CellSize };

	GDALRasterBand *poBand;

	poDstDS->SetGeoTransform(adfGeoTransform );

	poDstDS->SetProjection(SpatialRef.data());


	poBand = poDstDS->GetRasterBand(1);
	poBand->SetNoDataValue(m_fNullDataValue);

	poBand->RasterIO( GF_Write, 0, 0, NumberOfCols, NumberOfRows, m_pDataBuf, NumberOfCols, NumberOfRows, GDT_Float32, 0, 0 );    

	GDALClose( (GDALDatasetH) poDstDS );

}


float CRasterLayer::getXMin()
{
	return Xmin;
}
float CRasterLayer::getXMax()
{
	return Xmax;
}

float CRasterLayer::getYMin()
{
	return Ymin;
}

float CRasterLayer::getYMax()
{
	return Ymax;
}


float CRasterLayer::getNullDataValue()
{
	return m_fNullDataValue;
}


int CRasterLayer::rowOfAPoint(double y)
{ 
	
	return floor((float)NumberOfRows-(y-Ymin)/CellSize);
}

int CRasterLayer::colOfAPoint(double x)
{
	
	return floor((x-Xmin)/CellSize);

}

double CRasterLayer::xOfACell(int col)
{
	return Xmin+CellSize*(col+0.5);
}

double CRasterLayer::yOfACell(int row)
{
	return Ymin+CellSize*(NumberOfRows-row-0.5);
}

string CRasterLayer::getSpatialRef()
{
	return SpatialRef;
}

void CRasterLayer::setSpatialRef(string spatialRef)
{
	SpatialRef = spatialRef;
}

int CRasterLayer::getDataType()
{
	return m_iDataType;
}

void CRasterLayer::setDataType(int dataType)
{
	m_iDataType = dataType;
}


void CRasterLayer::setXMin(float value)
{
	Xmin = value;
}
void CRasterLayer::setXMax(float value)
{
	Xmax = value;
}

void CRasterLayer::setYMin(float value)
{
	Ymin = value;
}

void CRasterLayer::setYMax(float value)
{
	Ymax = value;
}

void CRasterLayer::setCellSize(float cellSize)
{
	CellSize = cellSize;
}


void CRasterLayer::setNumberOfColumns(int colNum)
{
	NumberOfCols = colNum;
}
void CRasterLayer::setNumberOfRows(int rowNum)
{
	NumberOfRows = rowNum;
}

void CRasterLayer::setNullDataValue(float value)
{
	m_fNullDataValue = value;
}