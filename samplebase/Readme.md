注意事项：

1. 本程序采用的推测土壤属性值的公式为：propertyVals[rowIdx][colIdx] = SumSimilarityValue/SumSimilarity;
也可以采用另一个公式propertyVals[rowIdx][colIdx] = MaxValue * MaxSimilarity + (1.0 - MaxSimilarity) * (SumSimilarityValue - MaxValue * MaxSimilarity) / (SumSimilarity - MaxSimilarity);
2. 使用本程序时，请保证所有环境因子的nodata数值都是一致的（-9999最好）
3. 本程序的运行需要配置openmpi并行计算环境和gdal环境（Linux）
	编译本程序的命令：mpic++ *.cpp -o samplebase -lgdal 其中samplebase名称随意
	运行本程序的命令：mpirun -n 进程数 程序路径 参数 
 
例如：

**version 1.1 2015-09-03 **

mpirun -n 4 samplebase "/home/admin/work/sampleWebservice/data2/geo.tif#/home/admin/work/sampleWebservice/data2/jiangshui.tif#/home/admin/work/sampleWebservice/data2/slope.tif#/home/admin/work/sampleWebservice/data2/dem.tif" "/home/admin/work/sampleWebservice/data2/training.csv" "Geology?Boolean#Climate?Gaussian#Terrain?Gaussian#Terrain?Gaussian" 0.5 "/home/admin/work/sampleWebservice/result/property" "/home/admin/work/sampleWebservice/result/uncertainty"

备注：Gaussian可换为Gower

**version 1.2 2016-01-25**


mpirun -n 4 samplebase "/home/admin/work/test/data/jiangshui.tif#/home/admin/work/test/data/geo.tif#/home/admin/work/test/data/dem.tif#/home/admin/work/test/data/slope.tif" "/home/admin/work/test/data/training.csv" "Climate?Gower#Geology?Boolean#Terrain?Gower#Terrain?Gower" 0.5 "/home/admin/work/test/result/propertyClay" "/home/admin/work/test/result/uncertaintyClay" "X" "Y" "Clay" "Limit" "Limit"

备注：X，Y，Clay分别是样点文件中坐标x,y和待推测属性的列名，综合方法可以换为Average

