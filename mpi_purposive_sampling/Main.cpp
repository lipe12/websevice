#include <iostream>
#include <mpi.h>
#include <string>
#include "mpi_FCM.h"
#include "RasterLayer.h"
#include "RasterDatabase.h"
#include "Utility.h"
#include <math.h>



using namespace std;

int main(int argc, char *argv[])
{

	int minClassNumber;
	int maxClassNumber;
	float error;
	char* inputfile;
	int maxIteration;
	float alpha;
	int numOfPnt;
	float distThreshold;
	char* resultDir;
    
	double startTime=MPI_Wtime();

	if(argc != 10)
	{  
		cout<<"[ERROR]: please input correct parameter!"<<endl;		
		return 0;


	}
	else
	{
		inputfile = argv[1];  	
		minClassNumber = atoi(argv[2]); 
		maxClassNumber = atoi(argv[3]);
		error = (float)atof(argv[4]);
		maxIteration = atoi(argv[5]);
		alpha = atof(argv[6]);
		numOfPnt = atoi(argv[7]);
		distThreshold = atof(argv[8]);
		resultDir = argv[9];  
	}

	if(minClassNumber >= maxClassNumber)
	{
		cout <<"[ERROR]: minimum cluster number should be smaller than maximum cluster number!" << endl;
		return 0;
	}
 
	double t1=MPI_Wtime();
	MPI_Init(&argc,&argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char buffer[1000];
	string str;

	

	cout << "rank" << rank << ": FCM starts!" << endl;
	str = resultDir;
	str =  str + "/fcm_result";
	makeDir(str);


	for(int i=minClassNumber; i<=maxClassNumber; i++)
	{
		char*  outputfile = new char[1000];	
			
		sprintf(buffer,"%s%s%d%s",resultDir,"/fcm_result/",i ,"_clusters");
		string str = buffer;	
		makeDir(str);
		sprintf(outputfile,"%s%s%d%s",resultDir, "/fcm_result/", i, "_clusters/Fuzzy.tif");	
		mpi_FCM  fcm(i,maxIteration,error); 
		fcm.FCM_Processing(inputfile,outputfile,error);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
	cout << "rank" << rank << ":FCM is finished!" << endl;
	double t2=MPI_Wtime();
	double computingTime = t2 - t1;

    if(rank != 0)
        exit(0);	

	cout << "Pattern extraction starts!"<<endl;

	int itemNum = maxClassNumber - minClassNumber + 1;

	str = resultDir;	
	str=  str + "/alpha_cut";
	makeDir(str);


	str = resultDir;	
	str=  str + "/patterns";
	makeDir(str);

	str = resultDir;	
	str=  str + "/patterns_all";
	makeDir(str);

	int nrow,ncol;
	float XMin, YMin, XMax, YMax;
	float CellSize;
	string spatialRef;
	float noDataValue;


	sprintf(buffer, "%s%s%d%s", resultDir, "/fcm_result/",  minClassNumber, "_clusters/FuzzyDegree_1.tif");
	string strMemberFile = buffer;

	CRasterLayer* tmpMemberFile = new CRasterLayer(strMemberFile);
	nrow = tmpMemberFile->getNumberOfRows();
	ncol = tmpMemberFile->getNumberOfColumns();
	CellSize = tmpMemberFile->getCellSize();
	XMin = tmpMemberFile->getXMin();
	YMin = tmpMemberFile->getYMin();
	XMax = tmpMemberFile->getXMax();
	YMax = tmpMemberFile->getYMax();
    spatialRef = tmpMemberFile->getSpatialRef();
	noDataValue = tmpMemberFile->getNullDataValue();

	delete tmpMemberFile;

	

	for(int i=minClassNumber; i<=maxClassNumber; i++)
	{
		float *alphaCutIntegrationBuf = new float[nrow*ncol];
		for(int irow=0;irow<nrow;irow++)
			for(int icol=0;icol<ncol;icol++)		 	 
			{		
				alphaCutIntegrationBuf[irow*ncol+icol] = 0;
			}
			for(int j=0; j<i; j++)
			{
				char buffer[1000];			  
				sprintf(buffer, "%s%s%d%s%d%s", resultDir, "/fcm_result/", i, "_clusters/FuzzyDegree_", j+1, ".tif");
				string strMemberFile = buffer;

				CRasterLayer * memberFile = new CRasterLayer(strMemberFile);

				int readflag = memberFile->readData();

				if(readflag > 0)
				{
					cout <<"[ERROR]: fail to read clusting result file." << endl;
					return 0;

				}
				else
				{	
					for(int irow=0;irow<nrow;irow++)
						for(int icol=0;icol<ncol;icol++)		 	 
						{

							if( memberFile->getData(irow, icol)>=alpha)
							{
								alphaCutIntegrationBuf[irow*ncol+icol] += 1*(j+1);
							}
							else 
								alphaCutIntegrationBuf[irow*ncol+icol] += 0*(j+1);			  
						}  

				}

				delete memberFile;
			}



			//set integrated alpha cut file
			sprintf(buffer, "%s%s%d%s", resultDir,"/alpha_cut/", i, "_alpha_cut.tif");		
			string strOutputFileName = buffer;
			CRasterLayer* output = new CRasterLayer();
			output->setFileName(strOutputFileName);
			output->setNumberOfRows(nrow);
			output->setNumberOfColumns(ncol);
			output->setXMin(XMin);
			output->setYMin(YMin);
			output->setXMax(XMax);
			output->setYMax(YMax);
			output->setSpatialRef(spatialRef);
			output->setCellSize(CellSize);
			output->setNullDataValue(noDataValue);
			output->m_pDataBuf = alphaCutIntegrationBuf;	
            
			output->wirteData("GTiff");
			delete output;			


	}



	CRasterDatabase *pAlphaCutDatabase = new CRasterDatabase();
	for(int i=minClassNumber; i<=maxClassNumber; i++)
	{

		sprintf(buffer, "%s%s%d%s",resultDir,"/alpha_cut/", i, "_alpha_cut.tif");	
		string strLayerPath = buffer;

		CRasterLayer *pLayer = new CRasterLayer(strLayerPath);
		pLayer->readData();
		pAlphaCutDatabase->AddLayer(pLayer);

	}




	int pattern[MAX_ITEM_NUM];	
	int stability[MAX_ITEM_NUM];
	int patternNum=0;
	for(int i=0;i<MAX_ITEM_NUM;i++)
	{
		pattern[i]=-3;  
		stability[i] = 0;
	}

	int num = 0;
	for(int irow=0;irow<nrow;irow++)
		for(int icol=0;icol<ncol;icol++)		 	 
		{	

			num++;

			int iDifferent=0;
			for(int n=0;n<patternNum;n++)
			{
				int matchNum = 0;
				int zeroNum = 0;
				int mismatchNum1 = 0;
				int mismatchNum2 = 0;
				int m;
				for(int i=0;i<itemNum;i++)
				{			

					if(fabs(pAlphaCutDatabase->getLayer(i)->getData(irow,icol) - 0) <VERYSMALL && pattern[n*itemNum+i] == 0 )					  
						zeroNum ++;
					else if(fabs(pAlphaCutDatabase->getLayer(i)->getData(irow,icol) - pattern[n*itemNum+i]) < VERYSMALL)
						matchNum ++;

					else if(fabs(pAlphaCutDatabase->getLayer(i)->getData(irow,icol) - 0 ) >= VERYSMALL && pattern[n*itemNum+i] == 0)				     
						mismatchNum1 ++;	

					else if(fabs(pAlphaCutDatabase->getLayer(i)->getData(irow,icol) - 0) < VERYSMALL && pattern[n*itemNum+i] != 0)
						mismatchNum2 ++;



				}

				if(zeroNum + matchNum + mismatchNum1 == itemNum && mismatchNum1 > 0 )  //this is higher level pattern
				{

					for(int i=0;i<itemNum;i++)	   
						pattern[n*itemNum+i] = int(pAlphaCutDatabase->getLayer(i)->getData(irow,icol));			  
					stability[n] = matchNum + mismatchNum1;

				}
				else if(zeroNum + matchNum + mismatchNum2 == itemNum && mismatchNum2 > 0 )  //this is lower level pattern
					m = 0; //do nothing


				else if(zeroNum + matchNum == itemNum)
					m = 0;  //do nothing
				else
					iDifferent++;
			}

			if(iDifferent == patternNum)  //this is totally different pattern
			{
				for(int i=0;i<itemNum;i++)	
				{
					pattern[patternNum*itemNum+i] = int(pAlphaCutDatabase->getLayer(i)->getData(irow,icol));	
					if(pattern[patternNum*itemNum+i] != 0)
						stability[patternNum]++; 
				}

				patternNum++;
			}
		}



		for(int i=0; i<patternNum; i++)
			for(int j=i+1; j<patternNum;j++)
			{
				bool tag = true;
				for(int m=0;m<itemNum;m++)	
				{
					if(pattern[i*itemNum+m] != pattern[j*itemNum+m])
					{ 
						tag = false;
						break;
					}

				}
				if(tag == true)  //the same pattern found, move
				{
					for(int k=j;k<patternNum-1;k++)
					{
						for(int m=0;m<itemNum;m++)	
						{
							pattern[k*itemNum+m] = pattern[(k+1)*itemNum+m];

						}
						stability[k] = stability[k+1]; 


					}
					patternNum--;
					j--;


				}
			}



			int area[MAX_ITEM_NUM];

			vector<string > vecAllDesignedPnts;
			vector<int > vecAllDesignedPntPattenID;


			for(int i=0; i<patternNum; i++)
			{
				
				if(i==6)
				{
					int a = 3;
				}


				float *resultBuf = new float[nrow*ncol];
				float *resultBufAll = new float[nrow*ncol];
				int *patchBuf = new int[nrow*ncol];


				area[i] = 0;
				CRasterDatabase *pFuzzyDatabase = new CRasterDatabase();
				//open fuzzy membership maps for this pattern
				for(int j=0; j<itemNum;j++)
				{
					if(pattern[i*itemNum+j]!=0)
					{
						sprintf(buffer, "%s%s%d%s%d%s", resultDir, "/fcm_result/", minClassNumber + j, "_clusters/FuzzyDegree_",pattern[i*itemNum+j], ".tif");
						string strLayerPath = buffer;

						CRasterLayer *pLayer = new CRasterLayer(strLayerPath);
						pLayer->readData();

						pFuzzyDatabase->AddLayer(pLayer);

					}

				}


				for(int irow=0; irow<nrow; irow++)
					for(int icol=0; icol<ncol; icol++)
					{	
						resultBufAll[irow*ncol+icol]=0;
						int temp = 0;
						for(int j=0; j<itemNum; j++)
						{  
							if(pattern[i*itemNum+j]!=0)
							{
								resultBufAll[irow*ncol+icol]+= pFuzzyDatabase->getLayer(temp)->getData(irow,icol);

								temp ++; 
							}
						}

						resultBufAll[irow*ncol+icol] = resultBufAll[irow*ncol+icol]/stability[i];


					}
					//output .tif for every pattern 	  
					sprintf(buffer,"%s%s%d%s", resultDir, "/patterns_all/pattern_all", i+1, ".tif");
					string strOutputFileName = buffer;
					CRasterLayer* outputAll = new CRasterLayer();
					outputAll->setFileName(strOutputFileName);
					outputAll ->setNumberOfRows(nrow);
					outputAll ->setNumberOfColumns(ncol);
					outputAll ->setXMin(XMin);
					outputAll ->setYMin(YMin);
					outputAll->setXMax(XMax);
					outputAll->setYMax(YMax);
					outputAll->setSpatialRef(spatialRef);
					outputAll ->setCellSize(CellSize);
					outputAll ->setNullDataValue(noDataValue);
					outputAll ->m_pDataBuf = resultBufAll;
					outputAll ->wirteData("GTiff");
					delete outputAll;		



					bool tag;
					//assign 
					for(int irow=0;irow<nrow;irow++)
						for(int icol=0;icol<ncol;icol++)	
						{
							tag = true;
							resultBuf[irow*ncol+icol]=0;

							int temp = 0;

							for(int j=0;j<itemNum;j++)	
							{

								if(pAlphaCutDatabase->getLayer(j)->getData(irow,icol)!=  pattern[i*itemNum+j])
								{
									tag = false;
									break;
								}
								else if(pattern[i*itemNum+j]!=0)
								{
									resultBuf[irow*ncol+icol]+= pFuzzyDatabase->getLayer(temp)->getData(irow,icol);

									temp ++; 
								}
							}
							if(tag == true)  //the pixel match with the pattern
							{
								area[i] ++;
								resultBuf[irow*ncol+icol] = resultBuf[irow*ncol+icol]/stability[i];

							}
							else 
							{
								resultBuf[irow*ncol+icol] = noDataValue;

							}
						}     



								float *highestValues = new float[numOfPnt];
								int *rowsHighest = new int[numOfPnt];
								int *colsHighest = new int[numOfPnt];



								for(int k=0; k<numOfPnt; k++)
								{
									float max = 0;
									int curRow = -1;
									int curCol = -1;


									for(int irow=0;irow<nrow;irow++)
										for(int icol=0;icol<ncol;icol++)
										{
											bool visited = false;

											for(int m=0; m<k; m++)
											{

												if(sqrt(pow(double(irow - rowsHighest[m]), 2) + pow(double(icol-colsHighest[m]), 2)) * CellSize < distThreshold)
												{
													visited = true;
													break;
												}
											}
											if(visited == false)
											{
												if(resultBuf[irow*ncol+icol] > max)
												{

													max = resultBuf[irow*ncol+icol];
													curRow = irow;
													curCol = icol;
												}
											}
										}



										highestValues[k]  = max;
										rowsHighest[k] = curRow;
										colsHighest[k] = curCol;




										if(curRow != -1)
										{

											string onePnt;
											sprintf(buffer, "%d,%d,%d,%f,%f,%f", stability[i], i+1 ,area[i], XMin + CellSize * (colsHighest[k]+0.5), YMin + CellSize * (nrow - rowsHighest[k]-0.5),  highestValues[k]);
											onePnt = buffer;
											vecAllDesignedPnts.push_back(onePnt);
											vecAllDesignedPntPattenID.push_back(i+1);
										}




								}


								sprintf(buffer, "%s%s%d%s", resultDir, "/patterns/pattern", i+1, ".tif");
								strOutputFileName = buffer;
								CRasterLayer* output = new CRasterLayer();
								output->setFileName(strOutputFileName);
								output->setNumberOfRows(nrow);
								output->setNumberOfColumns(ncol);
								output->setXMin(XMin);
								output->setYMin(YMin);
								output->setXMax(XMax);
								output->setYMax(YMax);
								output->setSpatialRef(spatialRef);
								output->setCellSize(CellSize);
								output->setNullDataValue(noDataValue);
								output->m_pDataBuf = resultBuf;
								output->wirteData("GTiff");
								delete output;			



								if(pFuzzyDatabase != NULL)
								{
									delete pFuzzyDatabase;
									pFuzzyDatabase = NULL;
								}

								if(patchBuf != NULL)
								{
									delete patchBuf;
									patchBuf = NULL;
								}


			}


			if(pAlphaCutDatabase != NULL)
			{
				delete pAlphaCutDatabase;
				pAlphaCutDatabase = NULL;
			}





			int patternIDs[1000];

			for(int i=0; i<patternNum; i++)
			{
				patternIDs[i] = i+1;
			}

			int r;

			int temp;

			for(int i=0; i<patternNum-1; i++)
			{
				r=i;
				for(int j=i+1; j<patternNum; j++)
				{
					if(stability[r]<stability[j])
						r = j;
					else if(stability[r]==stability[j] && area[r]<area[j])
						r = j;
				}
				if(r!=i)  //exchange
				{
					for(int m=0; m<itemNum; m++)
					{
						temp = pattern[i*itemNum+m];
						pattern[i*itemNum+m] = pattern[r*itemNum+m];
						pattern[r*itemNum+m]= temp;

					}
					temp=stability[i];
					stability[i] = stability[r];
					stability[r]= temp;

					temp=area[i];
					area[i] = area[r];
					area[r]= temp;


					temp=patternIDs[i];
					patternIDs[i] = patternIDs[r];
					patternIDs[r]= temp;


				}
			}


			cout<<"Pattern extraction is finished!"<<endl;

			cout<<"Write results to files."<<endl;

			FILE* fp;

			string reportFileName = resultDir;
			reportFileName = reportFileName + "/pattern_List.csv";
			

			//output csv file 
			if ((fp=fopen(reportFileName.c_str(), "w"))==NULL) 
			{
				cout << "[ERROR]: Can't open result file for writing." << endl;
				return 0;

			}

			fprintf(fp, "Stability,ID,Name,Area");


			for (int i=0; i<patternNum; i++)
			{
				string name = "";

				for(int m=0; m<itemNum; m++)
				{
					if(pattern[i*itemNum+m]!=0)
					{
						if(name!="")				 //this first
						{
							string str;
							sprintf(buffer, "%s%d%s%d", "-", minClassNumber + m, "Class", pattern[i*itemNum+m]);
							str = buffer;
							name += str;
						}
						else
						{
							string str;
							sprintf(buffer, "%d%s%d", minClassNumber + m, "Class", pattern[i*itemNum+m]);
							str = buffer;
							name += str;
						}

					}
				}


	          // string oneLine;
	          sprintf(buffer, "\n%d,%d,%s,%d",stability[i], patternIDs[i] ,name.c_str(),area[i]);
			  fprintf(fp, buffer);


			}

			fclose(fp);


			string pntFileName = resultDir;
			pntFileName = pntFileName+ "/designed_samples.csv";


			//output csv file 
			if ((fp=fopen(pntFileName.c_str(), "w"))==NULL) 
			{
				cout << "[ERROR]: Can't open sample file to write." << endl;
				return 0;

			}

			fprintf(fp, "Stability,PatternID,TotalArea,RecommendedX,RecommendedY,Ave.Membership");



			for (int i=0; i<patternNum; i++)
			{

				for(int k = 0; k < vecAllDesignedPntPattenID.size(); k++)
				{
					if (vecAllDesignedPntPattenID[k] == patternIDs[i])     
					{

						string str = "\n" + vecAllDesignedPnts[k];
						fprintf(fp,str.c_str());
					}

				}

			}


			fclose(fp);

			cout<<"Purposive sampling design is finished!"<<endl;

    double endTime=MPI_Wtime();
	double totalTime = endTime - startTime;


    cout<<"[DEBUG][TIMESPAN][IO]" << totalTime - computingTime << endl; 
    cout<<"[DEBUG][TIMESPAN][COMPUTING]" << computingTime << endl; 
    cout<<"[DEBUG][TIMESPAN][TOTAL]" << totalTime << endl; 

}
