/*算法描述:
  FCM聚类算法采用的是给定类的个数K,将N个元素(对象)分配到K个类中去使得类内对象之间的相似性最大,而类之间的相似性最小 */
#include "mpi_FCM.h"
#define INF_T 0.01
mpi_FCM::mpi_FCM(int number,int maxIteration,float error)
{
	this->number = number;
	this->maxIteration=maxIteration;
	//this->tolerance=0.65;     //2.5m
	//this->tolerance=0.45;  //1m
	this->tolerance=error;     //0.5m
	this->Fuzzyness=2.0;
	this->format = "GTiff";
	this->pTransform = new double[6];
	pMyarray = NULL;
	pClass = NULL;
	GDALAllRegister();
}

mpi_FCM::~mpi_FCM()
{
}

int mpi_FCM::FCM_Processing(char* inputfile,char* outputfile,float error)
{	
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs); //进程数
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);    //当前进程编号
	int i=0;
	int j=0;
	char strDest[2000];
	strcpy(strDest,inputfile);
	strParsing(strDest); //字符串解析	
	imageNum=allPathfiles.size(); //影像数目
	readInfo(allPathfiles[0]);
	double startTime=MPI_Wtime();
	for(i=0;i<imageNum;i++)
	{
		readDataFile(allPathfiles[i]);
		vecAllData.push_back(pMyarray);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	double endTime=MPI_Wtime();
	double IOReadTime=endTime-startTime;
	//cout<<"IO Read Time: "<<IOReadTime<<endl;
	
	center  =new double[number * imageNum];    //聚类中心 number * imageNum 
	int* centerIndex =new int [number];   //聚类中心索引 number * 1
	int* centerIndexTmp =new int [number];   //聚类中心索引 K * 1
	vector<int> tempIndex;
	for(i=0;i<subRow * column;i++)
	{
		if(vecAllData[0][i]!=noData && fabs(vecAllData[0][i]+9999)>INF_T)
			tempIndex.push_back(i);
	}
	CreateRandom(number,tempIndex.size(),centerIndexTmp);
	for(i=0;i<number;i++)
	{
		centerIndex[i]=tempIndex[centerIndexTmp[i]];
	}
    
	//2.5m迭代40次（单节点）
	if(rank==0)
	{
			for(i=0;i<number;i++)
		{
			for(j=0;j<imageNum;j++)
			{	
				center[i*imageNum+j]=vecAllData[j][centerIndex[i]];
			}
		}
		tempIndex.clear();
	}
	
	int curIteration=1;   //当前迭代次数
	sumNumerator=new double [number*imageNum];  //各进程中汇总各进程中分子求和
	sumDenominator=new double[number];  //各进程中汇总各进程中分母求和（隶属度求和）
	toalNumerator=new double [number*imageNum];  //主进程中汇总各进程中分子求和
	toalDenominator=new double[number];  //主进程中汇总各进程中分母求和（隶属度求和）
	degree=new double *[number];
	for(i=0;i<number;i++)
	{
		degree[i]=new double[block];
	}
	double subObjectValue;
	double subObjectValueNew;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(center, number*imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //主进程广播K个聚类中心
	InitDegree(); //初始化隶属度
	subObjectValueNew=fnObjectiveFunction();
	MPI_Allreduce(&subObjectValueNew,&objectValueNew,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	//objectValueNew=sqrt(objectValueNew);
	double df=0.0;
	MPI_Barrier(MPI_COMM_WORLD);
	do
	{

		objectValue=objectValueNew;
		fnClusterCentroid();  //计算聚类中心
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(sumNumerator,toalNumerator,number*imageNum,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(sumDenominator,toalDenominator,number,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		for(i=0;i<number;i++)
		{
			for(j=0;j<imageNum;j++)
			{
				center[i*imageNum+j]=toalNumerator[i*imageNum+j]/toalDenominator[i];
			}
		}

		fnDegreeMembership(); //计算隶属度
		
		subObjectValueNew=fnObjectiveFunction();
		
		MPI_Allreduce(&subObjectValueNew,&objectValueNew,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		
		//objectValueNew=sqrt(objectValueNew);
		df=fabs(objectValueNew-objectValue);
     		MPI_Barrier(MPI_COMM_WORLD);
		if(df<tolerance)
			break;
		curIteration++;
	}while(curIteration<=maxIteration);
	MPI_Barrier(MPI_COMM_WORLD);
	double allPartitionCoef=0.0;
	double allEntropy=0.0;
	MPI_Allreduce(&partitionCoef,&allPartitionCoef,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&entropy,&allEntropy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	startTime=MPI_Wtime();	
	
	char* token = strtok(outputfile,".");
	for(i=0;i<number+1;i++)
	{	
		char *pathfiles=new char[100];
		if(i==0)
		{
			sprintf(pathfiles,"%s.tif",token,i);
		}
		else
		{
			sprintf(pathfiles,"%sDegree_%d.tif",token,i);
		}
		outPathfiles.push_back(pathfiles);
	}



	if(rank==0)
	{
		char *centerpath=new char[100];
		sprintf(centerpath,"%s.txt",token,i);
		FILE *fp=fopen(centerpath,"w");
		fprintf(fp,"NumberOfIteration：%d\n",curIteration-1);
		fprintf(fp,"MaxError：%f\n",df);
		fprintf(fp,"PartitionCoefficient：%f\n",allPartitionCoef);
		fprintf(fp,"Entropy：%f\n",allEntropy);
		fprintf(fp,"PayOff：%f\n",objectValueNew);
		for(i=0;i<number;i++)
		{
			fprintf(fp,"class%d: ",i);
			for(j=0;j<imageNum;j++)
			{
				fprintf(fp,"%6.3f\t",center[i*imageNum+j]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		//cout<<"Iteration times is "<<curIteration<<endl;
		for(i=0;i<number+1;i++)
		{
		//	cout<<outPathfiles[i]<<endl;
			createFile(outPathfiles[i],i);
			
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
    int loopCode = 1;
	MPI_Status status;

    if (rank != 0)
        MPI_Recv(&loopCode, 1, MPI_INT,  rank-1, 0, MPI_COMM_WORLD, &status);

    //cout << "[DEBUG]\t" << rank << endl;
	for(i = 0; i < number+1; i++)
		writeDataFile(outPathfiles[i], i);
	
	if (rank != numProcs-1)
	    MPI_Send(&loopCode, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	endTime=MPI_Wtime();
	double IOWriteTime=endTime-startTime;
	//cout<<"IO Write Time: "<<IOWriteTime<<endl;

	delete []center;
	delete []centerIndex;	
	delete []centerIndexTmp;
	delete []sumNumerator;
	delete []sumDenominator;
	delete []toalNumerator;
	delete []toalDenominator;
	//delete []degree;
	delete []this->pTransform;
	for(i=0;i<number;i++)
	{
		delete []degree[i];
	}
	delete []degree;
	for(i=0;i<imageNum;i++)
	{
		free(vecAllData[i]);
	}
	vecAllData.clear();
	return 0;
}

void  mpi_FCM::strParsing(char* inputfile)   //字符串解析
{
	char* token = strtok(inputfile,",");
	char *pathfile;
	while(NULL != token)
	{
		pathfile=token;
		allPathfiles.push_back(pathfile);
		token=strtok(NULL,",");
	}
}

void mpi_FCM::readInfo(const char* filename) //解析信息
{
	GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	if( poDatasetsrc == NULL )
	{
		cout<<"[ERROR] data file is not open correct"<<endl;
		exit(1);
	}
	poDatasetsrc->GetGeoTransform(pTransform);
	projection = poDatasetsrc->GetProjectionRef();

	GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand( 1 );
	noData = poBandsrc->GetNoDataValue();

	row = poBandsrc->GetYSize();
	column = poBandsrc->GetXSize();

	interval = (row + numProcs - 1)/numProcs;
	double *pMyarrayTmp;
	if ((row-(numProcs-1)*interval) <= 0)
		interval -= 1;

	    if ( rank == (numProcs - 1) )
	   {
		subRow = row - (numProcs - 1)*interval;
		block = subRow*column;
		pMyarrayTmp = new double[block];
	//	pMyarrayTmp = (double*)CPLMalloc(sizeof(double)*block);
	   }  
	   else
	   {
		subRow = interval;
		block = subRow*column;
		pMyarrayTmp = new double[block];
		//pMyarrayTmp = (double*)CPLMalloc(sizeof(double)*block);
	   }

	if (pMyarrayTmp == NULL)
	{
		cout<<"[ERROR] the allocation of memory is error!"<<endl;
		MPI_Finalize();
		return;
	}

	poBandsrc->RasterIO(GF_Read, 0, rank*interval, column, subRow, pMyarrayTmp, column, subRow, GDT_Float64, 0, 0);

	if (poDatasetsrc != NULL)
	{
		GDALClose((GDALDatasetH)poDatasetsrc);
		poDatasetsrc = NULL;
	}
	
	int i=0;
	int j=0;
	int *bRowData=new int[subRow];
	int *rowData=new int[row];
	int subPixelNum=0;
	for(i=0;i<subRow;i++)
	{
		bRowData[i]=0;
		for(j=0;j<column;j++)
		{
			if(pMyarrayTmp[i*column+j]!=noData && fabs(pMyarrayTmp[i*column+j]+9999)>INF_T)
			{
				bRowData[i]++;
				subPixelNum++;
			}
		}	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&subPixelNum,&samplingPixelNum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	int *blockNum=new int[numProcs];
	MPI_Gather(&subRow,1,MPI_INT,blockNum,1,MPI_INT,0,MPI_COMM_WORLD);
	int *displs=new int[numProcs];
	displs[0]=0;
	for(i=1;i<numProcs;i++)
		displs[i]=displs[i-1]+blockNum[i-1];
	MPI_Gatherv(bRowData,subRow,MPI_INT,rowData,blockNum,displs,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(rowData, row, MPI_INT, 0, MPI_COMM_WORLD); 
	for(i=1;i<row;i++)
	{
		rowData[i]=rowData[i-1]+rowData[i];
	}
	subData=new int[numProcs];
	int *NumTmp=new int[numProcs];

	for(i=0;i<numProcs-1;i++)
	{
		NumTmp[i]=(rowData[row-1]+numProcs - 1)/numProcs *(i+1);
	}
	NumTmp[numProcs-1]=rowData[row-1];
	j=0;
	for(i=0;i<row-1;i++)
	{
		if(rowData[i]<NumTmp[j]&&rowData[i+1]>=NumTmp[j])
		{
			subData[j]=i+1;
			j++;
			
		}		
	}
	subData[numProcs-1]=row-1;
	MPI_Barrier(MPI_COMM_WORLD);
	int *subTmp=new int[numProcs];
	subTmp[0]=subData[0]+1;
	for(i=1;i<numProcs;i++)
	{
		subTmp[i]=subData[i]-subData[i-1];
	}
	subRow=subTmp[rank];
	for(i=numProcs-1;i>0;i--)
	{
		subData[i]=subData[i-1]+1;
	}
	subData[0]=0;
	delete []subTmp;
	delete []blockNum;
	delete []displs;

	delete []pMyarrayTmp;
	//free(pMyarrayTmp);
}

void mpi_FCM::readDataFile(const char* filename)
{
	GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	if( poDatasetsrc == NULL /*检查是否正常打开文件*/)
	{
		cout<<"[ERROR] data file is not open correct"<<endl;
		exit(1);
	}
	poDatasetsrc->GetGeoTransform(pTransform);
	projection = poDatasetsrc->GetProjectionRef();

	GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand( 1 );
	noData = poBandsrc->GetNoDataValue();
/*
	row = poBandsrc->GetYSize();
	column = poBandsrc->GetXSize();

	interval = (row + numProcs - 1)/numProcs;

	if ( rank == (numProcs - 1) )
	{
		subRow = row - (numProcs - 1)*interval;
		block = subRow*column;
		pMyarray = (double*)CPLMalloc(sizeof(double)*block);
		pClass = (int*)CPLMalloc(sizeof(int)*block);
	}
	else
	{
		subRow = interval;
		block = subRow*column;
		pMyarray = (double*)CPLMalloc(sizeof(double)*block);
		pClass = (int*)CPLMalloc(sizeof(int)*block);
	}

*/
	block = subRow*column;
	pMyarray = (double*)CPLMalloc(sizeof(double)*block);
	pClass = (int*)CPLMalloc(sizeof(int)*block);
	if (pMyarray == NULL)
	{
		cout<<"[ERROR] the allocation of memory is error!"<<endl;
		MPI_Finalize();
		return;
	}
	if (pClass == NULL)
	{
		cout<<"[ERROR] the allocation of memory is error!"<<endl;
		MPI_Finalize();
		return;
	}
	poBandsrc->RasterIO(GF_Read, 0, subData[rank], column, subRow, pMyarray, column, subRow, GDT_Float64, 0, 0);

	if (poDatasetsrc != NULL)
	{
		GDALClose((GDALDatasetH)poDatasetsrc);
		poDatasetsrc = NULL;
	}
}

void mpi_FCM::createFile(const char *destfilename,int tag)
{
	GDALDriver* poDriver = NULL;
	poDriver = GetGDALDriverManager()->GetDriverByName(format.c_str());
	if (poDriver == NULL)
	{
		cout<<"[ERROR] poDriver is NULL."<<endl;
		exit(1);
	}
	char **papszMetadata = poDriver->GetMetadata();
	GDALDataset* poDataset = NULL;
	if(tag==0)	
		poDataset = poDriver->Create(destfilename, column, row, 1, GDT_Int32, NULL);
	else
		poDataset = poDriver->Create(destfilename, column, row, 1, GDT_Float64, NULL);
	poDataset->SetGeoTransform( pTransform );	
	poDataset->SetProjection(projection.c_str());
	
	if (poDataset == NULL)
	{
		cout<<"[ERROR] poDatasetdest is NULL"<<endl;
		exit(1);
	}

	if (poDataset != NULL)
	{
		poDataset->FlushCache();
		GDALClose((GDALDatasetH)poDataset);
		poDataset = NULL;
	}
}

void mpi_FCM::writeDataFile(const char* filename,int tag)  //tag=0,输出结果，tag=1或其他，输出模糊隶属度
{
	GDALDataset* poDataset = NULL;
	poDataset = (GDALDataset *) GDALOpen( filename, GA_Update );
	if( poDataset == NULL /*检查是否正常打开文件*/)
	{
		cout<<"[ERROR] data file is not open correct"<<endl;
		exit(1);
	}
	GDALRasterBand*	poBanddest = poDataset->GetRasterBand(1);
	if (poBanddest == NULL)
	{
		cout<<"[ERROR] poBanddest is NULL"<<endl;
		exit(1);
	}
	if(rank == 0)
	{
		poBanddest->SetNoDataValue(this->noData);
	}
	if(tag==0)
	{
		poBanddest->RasterIO(GF_Write, 0,  subData[rank], column, subRow, pClass, column, subRow, GDT_Int32, 0, 0);
	}
	else
	{
		poBanddest->RasterIO(GF_Write, 0,  subData[rank], column, subRow, degree[tag-1], column, subRow, GDT_Float64, 0, 0);
	}
	poDataset->FlushCache();
	
	if (poDataset != NULL)
	{
		GDALClose((GDALDatasetH)poDataset);
		poDataset = NULL;
	}
}

void mpi_FCM::CreateRandom(int k,int n,int *centerIndex) //产生可以随输出控制的 k与n （可舍弃） 
{
	int i=0,j=0;
	srand((unsigned int)time(NULL));
	for(i=0;i<k;i++)
	{
		int a=rand()%n;
		for(j=0;j<i;j++)
		{
			if(centerIndex[j]==a)
				break;
		}
		if(j>=i)
		{
			centerIndex[i]=a;
		}
		else
		{
			i--;
		}
	}
} 

//初始化隶属度
void mpi_FCM::InitDegree()
{
	double diff;
	int i=0;
	int j=0;
    for (j = 0; j < number; j++)
    {
        for (i = 0;i < block;i++)
        {
			if(vecAllData[0][i]==noData || fabs(vecAllData[0][i]+9999)<=INF_T)
			{
				degree[j][i]=noData;
				pClass[i]=noData;
			}
			else
			{
				diff = fnDistance(i,j);
				//degree[j][i] = (diff == 0.0) ? Eps:diff;
				degree[j][i] = (diff == 0.0) ? 1.0:1.0/diff;
			}
        }
    }
    fnClusterMembershipValues();
}

 //使每一个样本与各个聚类中心的隶属度之和为1
void mpi_FCM::fnClusterMembershipValues()
{
	int i=0;
	int j=0;
    for (i = 0; i < block; i++)
    {
		if(vecAllData[0][i]!=noData && fabs(vecAllData[0][i]+9999)>INF_T)
		{
			double max = 0.0;
			double min = pow((float)10,(float)6);
			double sum = 0.0;
			double newmax = 0.0;
			pClass[i]=0;
			for (j = 0; j < number; j++)
			{
				max = degree[j][i] > max ? degree[j][i]:max;
				min = degree[j][i] < min ? degree[j][i]:min;
			}
			//使样本与聚类中心的隶属度在 0 和 1之间
			for (j = 0; j < number; j++)
			{
				if(max==min)
				{
					degree[j][i]=1.0/number;
				}
				else
				{
					degree[j][i] = (degree[j][i] - min) / (max - min);	
				}
				sum += degree[j][i];
			}
			//使每一个样本与各个聚类中心的隶属度之和为1
			for (j = 0; j < number; j++)
			{	
				if(sum!=0.0)
				  degree[j][i] = degree[j][i] / sum;
				/*if (double.IsNaN(degree[i][j]))
				{
					degree[i][j] = 0.0;
				}*/
				if(degree[j][i]>newmax)
				{
					newmax=degree[j][i];
					pClass[i]=j;
				}
			   // newmax = degree[i][j] > newmax ? degree[i][j] : newmax;
			}
		}
    }
}

//计算距离
double mpi_FCM::fnDistance(int i, int j)  //第i个点与第j个中心的距离
{
	int img=0;
	double distance=0.0;
	for(img=0;img<imageNum;img++)
	{
		distance+=(vecAllData[img][i]-center[j*imageNum+img])*(vecAllData[img][i]-center[j*imageNum+img]);
	}
    return sqrt(distance);
}

//计算隶属度
void mpi_FCM::fnDegreeMembership()
{
	int i=0;
	int j=0;
	int k=0;
    for (j = 0; j < number; j++)
    {
        for (i = 0; i < block;i++)
        {
			if(vecAllData[0][i]!=noData && fabs(vecAllData[0][i]+9999)>INF_T)
			{
				double temp = 0.0;
				temp = fnDistance(i,j);
				if(temp == 0.0)
				{
					temp = Eps;
				}
				double sumDistance = 0.0;  //每一个样本到各个聚类中心的距离之和
				for (k = 0; k < number;k++ )
				{
					if(fnDistance(i,k)==0)
						sumDistance += pow(temp / Eps,2/(Fuzzyness-1));
					else
						sumDistance += pow(temp / fnDistance(i,k),2/(Fuzzyness-1));
				}
				degree[j][i] = 1.0 / sumDistance;
			}
        }
    }
    fnClusterMembershipValues();
}

//计算聚类中心
void mpi_FCM::fnClusterCentroid()
{
	int i=0;
	int j=0;
	int img=0;
    for (j=0;j<number;j++)
    {
        double m = 0.0;
		for(img=0;img<imageNum;img++)
		{
			sumNumerator[j*imageNum+img]=0.0;  //分子求和
		}
        sumDenominator[j] = 0.0;   //分母求和
        for (i=0;i<block;i++)
        {
			if(vecAllData[0][i]!=noData && fabs(vecAllData[0][i]+9999)>INF_T)
			{
				double df = pow(degree[j][i],Fuzzyness);
				for(img=0;img<imageNum;img++)
				{
					sumNumerator[j*imageNum+img]+=df*vecAllData[img][i];
				}
				sumDenominator[j] += df; 
			}
        }
    }
}

//计算目标函数
double mpi_FCM::fnObjectiveFunction()
{
    double Jk = 0.0;
	int i=0;
	int j=0;
	int img=0;
    for (i = 0; i < block; i++)
    {
		if(vecAllData[0][i]!=noData && fabs(vecAllData[0][i]+9999)>INF_T)
		{
			for (j = 0; j < number; j++)
			{
				Jk += pow(degree[j][i], Fuzzyness) * fnDistance(i,j)*fnDistance(i,j);
			}
		}
    }
	
	//计算 熵和划分系数
	partitionCoef = 0.0;
	entropy = 0.0;
	for (i = 0; i < block; i++)
    {
		if(vecAllData[0][i]!=noData && fabs(vecAllData[0][i]+9999)>INF_T)
		{
			for (j = 0; j < number; j++)
			{
				partitionCoef += degree[j][i]*degree[j][i]/samplingPixelNum;
				if(degree[j][i]>0) 
				{
					entropy += -1.0*degree[j][i]*log(degree[j][i])/samplingPixelNum;
				}
			}
		}
    }	
    return Jk;
}
