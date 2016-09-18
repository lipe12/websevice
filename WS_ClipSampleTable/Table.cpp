#include "Table.h"
#include<fstream>

using namespace std;

Table::Table()
{
	Table::ColCount = 0;
	Table::RowCount = 0;
	return;
}

Table::Table(string filename)
{
	RowCount = 0;
	ColCount = 0;
	
	try
	{
		FILE *pfin = NULL;
		pfin = fopen(filename.c_str(),"r");
		if(pfin == NULL)
		{
			cerr<<"fail to open sample file"<<endl;
		}
		char row[500];
		if (fgets(row,500,pfin) != NULL)
		{
			string line = string(row);
			parseStr(line, ',', this->Header);
			ColCount = this->Header.size();
		}
		for (int i = 0; i < ColCount; i++)
		{
			this->_ColNames.push_back(this->Header[i]);
		}
		vector<string> tmprow;
		while (fgets(row,500,pfin) != NULL)
		{
			string line = string(row);
			parseStr(line, ',', tmprow);
			this->Values.push_back(tmprow);
			RowCount++;
			tmprow.clear();
		}

	}
	catch(double)
	{

	}

}

void Table::parseStr(string str, char c, vector<string> & tokens)
{
	unsigned int posl = 0;
	unsigned int posR = 0;
	while (posR < str.length()-1)
	{
		posR = str.find_first_of(c, posl);
		string sub = str.substr(posl, posR - posl);
		tokens.push_back(sub);
		posl = posR + 1;
	}
}

Table::Table(vector<string> Header, vector<vector<string> > Values)
{
	Table::RowCount = 0;
	Table::ColCount = 0;
	Table::Header = Header;
	for (unsigned int i = 0; i < Table::Header.size(); i++)
	{
		this->_ColNames.push_back(this->Header[i]);
		this->ColCount = this->ColCount + 1;
	}
	this->Values = Values;
	for (unsigned int i = 0; i < this->Values.size(); i++)
	{
		this->RowCount = this->RowCount + 1;
	}
}

vector<string> Table::GetHeader()
{
	return this->Header;
}

vector<vector<string> > Table::GetValues()
{
	return this->Values;
}

int Table::GetColumnIndexByName(string colName)
{
	int colIdx = -1;
	for(unsigned int i = 0; i < this->_ColNames.size(); i++)
	{
		if(colName == this->_ColNames[i])
			colIdx = i;
	}
	return colIdx;
}


string Table::GetValue(int Col, int Row)
{
	return this->Values[Row][Col];
}

void Table::SetValue(int Col, int Row, string Value)
{
	Values[Row][Col] = Value;
}

string Table::GetColName(int ColIndex)
{
	return this->_ColNames[ColIndex];
}

void Table::WriteAsCSVFile(string FileName)
{
	/*FILE *fp = fopen(FileName.c_str(), "W+");
	string HeaderString;
	HeaderString = this->Header[0];
	for (unsigned int i = 0; i < this->Header.size(); i++)
	{
		HeaderString = HeaderString + "," + this->Header[i];
	}
	fputs(HeaderString.c_str(), fp);
	for(unsigned int i = 0; i < this->Values.size(); i++)
	{
		string row = this->Values[i][0];
		for (unsigned int j = 1; j < this->Values[i].size(); j++)
		{
			row  = row + "," + this->Values[i][j];
		}
		fputs(row.c_str(), fp);
	}
	fclose(fp);*/


	ofstream fsp(FileName.c_str());
	if (fsp)
	{
		string HeaderString;
	    HeaderString = this->Header[0];
		fsp << HeaderString;
	    for (unsigned int i = 1; i < this->Header.size(); i++)
	   {
		   fsp<< "," << this->Header[i] ;
	    }
	    //fsp << endl;
	    for(unsigned int i = 0; i < this->Values.size(); i++)
	   {
		  // string row = this->Values[i][0];
		   fsp << Values[i][0];
		   for (unsigned int j = 1; j < this->Values[i].size(); j++)
		   {
			   fsp <<  "," << this->Values[i][j];
		    }
		   //fsp << endl;
	   }
		
		fsp.close();
	}
}