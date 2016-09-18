#pragma once

#include <string>
#include "RasterLayer.h"
using std::string;


//this class record basic setting of env database
class CSetting
{

public:
	CSetting(CRasterLayer *theData);
	CSetting(CSetting* setting);
	~CSetting();

private:
	int m_iColNum;
	int m_iRowNum;
	double m_dWest;
	double m_dEast;
	double m_dSouth;
	double m_dNorth;
	float m_fCellSize; 
	string m_strGridUnit;
	string m_strDataUnit;
	float m_fNodata;
    double m_dWestMin;
	double m_dEastMin;
	double m_dSouthMin;
	double m_dNorthMin;

	string m_strSpatialRef;


public:
	int getColNum();
	void setColNum(int iColNum);
	int getRowNum();
	void setRowNum(int iRowNum);
	double getWest();
	void setWest(double dWest);
	double getEast();
	void setEast(double dEast);
	double getSouth();
	void setSouth(double dSouth);
	double getNorth();
	void setNorth(double dNorth);
	double getCellSize();
	void setCellSize(float fCellSize);
	/*string getGridUnit();
	void setGridUnit(string strGridUnit);
	string getDataUnit();
	void setDataUnit(string strDataUnit);*/
	float getNodata();
	void setNodata(float fNodataMin);
	double getWestMin();
	void setWestMin(double dWestMin);
	double getEastMin();
	void setEastMin(double dEastMin);
	double getSouthMin();
	void setSouthMin(double dSouthMin);
	double getNorthMin();
	void setNorthMin(double dNorthMin);

	string getSpatialRef();
	void setSpatialRef(string spatialRef);




};