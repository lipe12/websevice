#ifndef TABLE_H
#define TABLE_H
#include<iostream>
#include<string>
#include<fstream>
#include<list>
#include<vector>
using namespace std;
class Table
{
public:
	//variables
	int ColCount;
	int RowCount;
	vector<string> Header;
	vector<vector<string> > Values;
	vector<string> _ColNames;
	//inits
	Table();
	Table(string filenames);
	Table(vector<string> Header, vector<vector<string> > Values);
	//Method
	vector<string> GetHeader();
	vector<vector<string> > GetValues();
	int GetColumnIndexByName(string colName);
	string GetValue(int Col, int Row);
	void SetValue(int Col, int Row, string Value);
	string GetColName(int ColIndex);
	void WriteAsCSVFile(string Filename);
	void parseStr(string str, char c, vector<string> & tokens);
};
#endif