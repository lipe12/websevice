#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "RasterLayer.h"
#include "Setting.h"
using namespace std;
#include <vector>

//class for environmental database
class CRasterDatabase
{
public:
	CRasterDatabase();
	CRasterDatabase(CRasterDatabase & anRasterDatabase);
	~CRasterDatabase();

public:
	CSetting* m_pTheSetting;  //setting records the extent, row,colun of the database

private:
	
	int m_iLayerNum;   //layer number
	vector <CRasterLayer *> m_vecEnvLayer; //vector to record data layers
	int m_iDEMLayerID; //record DEM layer ID in order to calculate surface distance
	bool m_bIsEditVersion;  //if the database is the edit version 
	CRasterDatabase *m_pBaseRasterDatabase; //if the database is edit version, this varibable records its base database pointer

public:
	int getLayerNum();
	void setLayerNum(int iLayerNum);
	int getDEMLayerID();
	CRasterLayer * getLayer(int iIndex);
	bool AddLayer(CRasterLayer* pLayer);
	void DeleteLayer(int iIndex);
	void DeleteAll();
	CRasterLayer* FindDataLayer(string strLayerName);  //find layer by name
	int FindDataLayerID(string strLayerName); //get the the order of layer in layer vector by matching layer name	

	void swapLayer(int iIndex1, int iIndex2);

};

