#pragma once

#include <string>
using std::string;


class CRasterLayer
{
public:
	CRasterLayer();	
	CRasterLayer(CRasterLayer& envLayer);
	CRasterLayer(string strLayerName, string strFileName);
	CRasterLayer(string strFileName);
	~CRasterLayer();

public:
	float *m_pDataBuf;

private:

	string m_strLayerName;  //layer name
	string m_strFileName;  //layer file
	int m_iDataType; //1: ratio/interval 2: norminal/categorical

	int NumberOfRows;
	int NumberOfCols;
	float CellSize;
	float Xmin;
	float Ymin;
	float Xmax;
	float Ymax;
	float m_fNullDataValue;
	string SpatialRef;

public:

	string getLayerName();
	void setLayerName(string strLayerName);

	string getFileName();
	void setFileName(string strFileName);	

	int getNumberOfColumns();
	int getNumberOfRows();

	void setNumberOfColumns(int colNum);
	void setNumberOfRows(int rowNum);


	string getSpatialRef();
	void setSpatialRef(string spatialRef);


	float getData(int iRow, int iCol);	
	float getData(double x,double y);

	int readHeader();
	int readData(); //0: success 1: fail
	int wirteData(string format); //0: success 1: fail

	int releaseDataBuf();

	float getCellSize();
    void setCellSize(float cellSize);


	float getXMin();
	float getXMax();
	float getYMin();
	float getYMax();


	void setXMin(float value);
	void setXMax(float value);
	void setYMin(float value);
	void setYMax(float value);
	
	float getNullDataValue();
	void setNullDataValue(float value);
	
	int rowOfAPoint(double y);    
	int colOfAPoint(double x);    
	double xOfACell(int col);
	double yOfACell(int row);


	int getDataType();
	void setDataType(int dataType);

};
