#include "RasterDatabase.h"
#include <iostream>
using std::string;


CRasterDatabase::CRasterDatabase()
{
	m_iLayerNum = 0;
	m_pTheSetting = NULL;
	m_iDEMLayerID = -1;
	m_bIsEditVersion = false;
	m_pBaseRasterDatabase = NULL;
}

//copy construction function
CRasterDatabase::CRasterDatabase(CRasterDatabase & anRasterDatabase)
{
	m_iLayerNum = 0;

	if(anRasterDatabase.m_pTheSetting!=NULL)
		m_pTheSetting = new CSetting(*(anRasterDatabase.m_pTheSetting));
	else
		m_pTheSetting = NULL;	

	for(int i=0;i<anRasterDatabase.m_iLayerNum;i++)
	{
		CRasterLayer *pTempLayer = new CRasterLayer(*(anRasterDatabase.getLayer(i)));
		AddLayer(pTempLayer);
	}
}



CRasterDatabase::~CRasterDatabase()
{
	if(m_pTheSetting != NULL)
	{
		delete m_pTheSetting;
		m_pTheSetting = NULL;
	}

	if(m_bIsEditVersion==true)  //if this is edit version, only remove layers, do not destroy them
	{
		while(m_iLayerNum != 0)
		{
			m_vecEnvLayer.pop_back();
			m_iLayerNum--;
		}
	}

	if(m_bIsEditVersion==false)  //if this is base version, destroy layers
	{
		while(m_iLayerNum != 0)
		{
			CRasterLayer *pTempLayer = m_vecEnvLayer[m_iLayerNum-1];
			m_vecEnvLayer.pop_back();
			if(pTempLayer != NULL)
			{
				delete pTempLayer;
				pTempLayer = NULL;
			}
			m_iLayerNum--;
		}
	}    
}

bool CRasterDatabase::AddLayer(CRasterLayer* pLayer)
{
	if(pLayer->getLayerName() == "elevation" ||pLayer->getLayerName() == "Elevation" ||pLayer->getLayerName() == "ELEVATION")
		m_iDEMLayerID = m_iLayerNum;

	if(m_iLayerNum == 0)
	{
		m_pTheSetting = new CSetting(pLayer);
		//m_pTheSetting->setSpatialRef(pLayer->getSpatialRef());

	}
	else
	{
		if(m_pTheSetting->getSpatialRef() != pLayer->getSpatialRef())
		{
			cout << "The spatial reference of the new layer does not match with the previous layers. The layer will not be added." <<endl;
			return false;

		}
		double prevWest,prevEast,prevSouth,prevNorth;	
		prevWest = m_pTheSetting->getWestMin();
		prevEast = m_pTheSetting->getEastMin();
		prevSouth = m_pTheSetting->getSouthMin();
		prevNorth = m_pTheSetting->getNorthMin();

		if(m_pTheSetting->getWestMin()< pLayer->getXMin())			
			m_pTheSetting->setWestMin(pLayer->getXMin());

		if(m_pTheSetting->getEastMin() > pLayer->getXMax())		
			m_pTheSetting->setEastMin(pLayer->getXMax());

		if(m_pTheSetting->getSouthMin() < pLayer->getYMin())			
			m_pTheSetting->setSouthMin(pLayer->getYMin());

		if(m_pTheSetting->getNorthMin() > pLayer->getYMax())
			m_pTheSetting->setNorthMin(pLayer->getYMax());

		if(m_pTheSetting->getWestMin() >= m_pTheSetting->getEastMin() || m_pTheSetting->getSouthMin() > m_pTheSetting->getNorthMin())
		{
			//return to previous state
			m_pTheSetting->setWestMin(prevWest);
			m_pTheSetting->setEastMin(prevEast);
			m_pTheSetting->setSouthMin(prevSouth);
			m_pTheSetting->setNorthMin(prevNorth);
			cout << "The spatial extent of the new layer does not have overlap with the previous layers. The layer will not be added." <<endl;
			return false;
		}


	}
	//pLayer->m_pTheData->calcStat();
	if(m_bIsEditVersion==true)
	{
		m_pBaseRasterDatabase->m_vecEnvLayer.push_back(pLayer);       		
		m_pBaseRasterDatabase->m_iLayerNum++; 
		m_vecEnvLayer.push_back(pLayer);
		m_iLayerNum++; 
	}
	else
	{
		m_vecEnvLayer.push_back(pLayer);
		m_iLayerNum++; 
	}

	return true;
}


void CRasterDatabase::DeleteLayer(int iIndex)
{
	//if the deleleted layer is the first layer, reset the setting
	if(iIndex ==0&&m_iLayerNum>1)
	{
		delete m_pTheSetting;
		m_pTheSetting = NULL;		
		CRasterLayer *pFirstLayer = m_vecEnvLayer[1];
		m_pTheSetting = new CSetting(pFirstLayer);
	}


	if(m_bIsEditVersion==true)
	{
		m_vecEnvLayer.erase(m_vecEnvLayer.begin()+iIndex);
		m_iLayerNum--;
	}

	else
	{
		CRasterLayer *pTempLayer = m_vecEnvLayer[iIndex];
		if(pTempLayer != NULL)
		{
			delete pTempLayer;
			pTempLayer = NULL;
		}
		m_vecEnvLayer.erase(m_vecEnvLayer.begin()+iIndex);
		m_iLayerNum--;

	}

	if(m_iLayerNum == 0)
	{
		delete m_pTheSetting;
		m_pTheSetting = NULL;
		return;
	}


	//recalculate the bounding box
	m_pTheSetting->setWest(m_vecEnvLayer[0]->getXMin());
	m_pTheSetting->setEast(m_vecEnvLayer[0]->getXMax());
	m_pTheSetting->setSouth(m_vecEnvLayer[0]->getYMin());
	m_pTheSetting->setNorth(m_vecEnvLayer[0]->getYMax());

	m_pTheSetting->setColNum(m_vecEnvLayer[0]->getNumberOfColumns());
	m_pTheSetting->setRowNum(m_vecEnvLayer[0]->getNumberOfRows());
	m_pTheSetting->setCellSize(m_vecEnvLayer[0]->getCellSize());

	//string gridUnit,dataUnit;
	//gridUnit.Format("%s",m_vecEnvLayer[0]->getGridUnits());
	//dataUnit.Format("%s",m_vecEnvLayer[0]->getDataUnits());
	//m_pTheSetting->setGridUnit(gridUnit);
	//m_pTheSetting->setDataUnit(dataUnit);

	m_pTheSetting->setNodata(m_vecEnvLayer[0]->getNullDataValue());

	m_pTheSetting->setWestMin(m_vecEnvLayer[0]->getXMin());
	m_pTheSetting->setEastMin(m_vecEnvLayer[0]->getXMax());
	m_pTheSetting->setSouthMin(m_vecEnvLayer[0]->getYMin());
	m_pTheSetting->setNorthMin(m_vecEnvLayer[0]->getYMax());


	for(int i=1; i<m_iLayerNum; i++)
	{
		if(m_pTheSetting->getWestMin()< m_vecEnvLayer[i]->getXMin())
			m_pTheSetting->setWestMin( m_vecEnvLayer[i]->getXMin());
		if(m_pTheSetting->getEastMin() > m_vecEnvLayer[i]->getXMax())
			m_pTheSetting->setEastMin( m_vecEnvLayer[i]->getXMax());
		if(m_pTheSetting->getSouthMin() < m_vecEnvLayer[i]->getYMin())
			m_pTheSetting->setSouthMin(m_vecEnvLayer[i]->getYMin());
		if(m_pTheSetting->getNorthMin() > m_vecEnvLayer[i]->getYMax())
			m_pTheSetting->setNorthMin(m_vecEnvLayer[i]->getYMax());

	}

}

void CRasterDatabase::DeleteAll()
{
	if(m_bIsEditVersion==true)
	{
		while(m_iLayerNum != 0)
		{
			m_vecEnvLayer.pop_back();
			m_iLayerNum--;
		}
	}
	else
	{
		while(m_iLayerNum != 0)
		{
			CRasterLayer *pTempLayer = m_vecEnvLayer[m_iLayerNum-1];
			m_vecEnvLayer.pop_back();
			if(pTempLayer != NULL)
			{
				delete pTempLayer;
				pTempLayer = NULL;
			}
			m_iLayerNum--;
		}
	}

}


int CRasterDatabase::getLayerNum()
{
	return m_iLayerNum;
}

void CRasterDatabase::setLayerNum(int iLayerNum)
{
	m_iLayerNum = iLayerNum;
}


CRasterLayer* CRasterDatabase::getLayer(int iIndex)
{
	if(iIndex>=0&&iIndex<m_iLayerNum)
		return m_vecEnvLayer[iIndex];
	else
		return NULL;
}

int CRasterDatabase::getDEMLayerID()
{
	return m_iDEMLayerID;
}


CRasterLayer*  CRasterDatabase::FindDataLayer(string strLayerName)
{
	int iLayerNum = getLayerNum();
	int i = 0;
	for(i=0; i<iLayerNum; i++)
	{
		if(getLayer(i)->getLayerName() == strLayerName)
		{
			return getLayer(i);
			break;
		}
	}
	if(i==iLayerNum)
		return NULL;

	return NULL;

}

int CRasterDatabase::FindDataLayerID(string strLayerName)
{
	int iLayerNum = getLayerNum();
	int i = 0;
	for(i=0; i<iLayerNum; i++)
	{
		if(getLayer(i)->getLayerName() == strLayerName)
		{
			return i;
			break;
		}
	}
	if(i==iLayerNum)
		return -1;

	return -1;

}


void CRasterDatabase::swapLayer(int iIndex1, int iIndex2)
{
	if(iIndex1>=this->getLayerNum()||iIndex2>=this->getLayerNum())
		return;
	CRasterLayer *layer1 = this->getLayer(iIndex1);
	CRasterLayer *layer2 = this->getLayer(iIndex2);

	this->m_vecEnvLayer.erase(m_vecEnvLayer.begin()+iIndex1);
	this->m_vecEnvLayer.insert(m_vecEnvLayer.begin()+iIndex1, layer2);
	this->m_vecEnvLayer.erase(m_vecEnvLayer.begin()+iIndex2);
	this->m_vecEnvLayer.insert(m_vecEnvLayer.begin()+iIndex2, layer1);



}


