#include "Setting.h"
#include "RasterLayer.h"

//The Setting class
CSetting::CSetting(CRasterLayer *theData)
{
	//CRasterLayer* theData = new CRasterLayer(aDataFile);
	
	m_iColNum = theData->getNumberOfColumns();
	m_iRowNum = theData->getNumberOfRows();
	m_fCellSize = theData->getCellSize(); 
	m_dWest = theData->getXMin();
	m_dEast = m_dWest+m_fCellSize*m_iColNum;
	m_dSouth = theData->getYMin();
	m_dNorth = m_dSouth+m_fCellSize*m_iRowNum;
	m_dWestMin = theData->getXMin();
	m_dEastMin = m_dWest+m_fCellSize*m_iColNum;
	m_dSouthMin = theData->getYMin();
	m_dNorthMin = m_dSouth+m_fCellSize*m_iRowNum;
	m_strSpatialRef = theData->getSpatialRef();
	//m_strGridUnit.Format("%s",theData->getGridUnits());
	//m_strDataUnit.Format("%s",theData->getDataUnits());
	m_fNodata = theData->getNullDataValue();
	//delete theData;
}



CSetting::CSetting(CSetting *aSetting)
{
	m_iColNum = aSetting->getColNum();
	m_iRowNum = aSetting->getRowNum();
	m_fCellSize = aSetting->getCellSize(); 
	m_dWest = aSetting->getWest();
	m_dEast = aSetting->getEast();
	m_dNorth = aSetting->getNorth();
	m_dSouth = aSetting->getSouth();
	m_dWestMin = aSetting->getWestMin();
	m_dEastMin = aSetting->getEastMin();
	m_dNorthMin = aSetting->getNorthMin();
	m_dSouthMin = aSetting->getSouthMin();
	//m_strGridUnit = aSetting->getGridUnit();
	//m_strDataUnit = aSetting->getDataUnit();	
	m_fNodata = aSetting->getNodata();
	m_strSpatialRef = aSetting->getSpatialRef();
}
CSetting::~CSetting()
{
}

int CSetting::getColNum()
{
	return m_iColNum;
}
void CSetting::setColNum(int iColNum)
{
	m_iColNum = iColNum;
}
int CSetting::getRowNum()
{
	return m_iRowNum;
}
void CSetting::setRowNum(int iRowNum)
{
	m_iRowNum = iRowNum;
}
double CSetting::getWest()
{
	return m_dWest;
}
void CSetting::setWest(double dWest)
{
	m_dWest = dWest;
}
double CSetting::getEast()
{
	return m_dEast;
}
void CSetting::setEast(double dEast)
{
	m_dEast = dEast;
}
double CSetting::getSouth()
{
	return m_dSouth;
}
void CSetting::setSouth(double dSouth)
{
	m_dSouth = dSouth;
}
double CSetting::getNorth()
{
	return m_dNorth;
}
void CSetting::setNorth(double dNorth)
{
	m_dNorth = dNorth;
}
double CSetting::getCellSize()
{
	return m_fCellSize;
}
void CSetting::setCellSize(float fCellSize)
{
	m_fCellSize = fCellSize;
}

float CSetting::getNodata()
{
	return m_fNodata;
}
void CSetting::setNodata(float fNodata)
{
	m_fNodata = fNodata;
}
double CSetting::getWestMin()
{
	return m_dWestMin;
}
void CSetting::setWestMin(double dWestMin)
{
	m_dWestMin = dWestMin;
}
double CSetting::getEastMin()
{
	return m_dEastMin;
}
void CSetting::setEastMin(double dEastMin)
{
	m_dEastMin = dEastMin;
}
double CSetting::getSouthMin()
{
	return m_dSouthMin;
}
void CSetting::setSouthMin(double dSouthMin)
{
	m_dSouthMin = dSouthMin;
}
double CSetting::getNorthMin()
{
	return m_dNorthMin;
}
void CSetting::setNorthMin(double dNorthMin)
{
	m_dNorthMin = dNorthMin;
}

string CSetting::getSpatialRef()
{
	return m_strSpatialRef;

}

void CSetting::setSpatialRef(string spatialRef)
{
	m_strSpatialRef = spatialRef;

}