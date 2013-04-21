MetaClusterTA:  MetaClusterTA.cpp Methods.o MetaCluster.o TomAlgorithm.o BaseStr.o ULLN.o KMER.o ShortKMER.o USet.o Reader.o Utils.o GLOBAL.o BWTs.o BWTDtStr.o 
	g++ -O2 -fopenmp MetaClusterTA.cpp Methods.o MetaCluster.o TomAlgorithm.o BaseStr.o ULLN.o KMER.o ShortKMER.o USet.o Reader.o Utils.o GLOBAL.o BWTs.o BWTDtStr.o -o MetaClusterTA

Methods.o: Methods.h Methods.cpp
	g++ -O2 -fopenmp Methods.h Methods.cpp -c

MetaCluster.o: MetaCluster.cpp MetaCluster.h
	g++ -O2 -fopenmp MetaCluster.cpp MetaCluster.h -c

Reader.o: Reader.h Reader.cpp
	g++ -O2 -fopenmp Reader.h Reader.cpp -c

TomAlgorithm.o: TomAlgorithm.cpp TomAlgorithm.h
	g++ -O2 -fopenmp TomAlgorithm.cpp  TomAlgorithm.h -c

ShortKMER.o: ShortKMER.h ShortKMER.cpp
	g++ -O2 -fopenmp ShortKMER.h ShortKMER.cpp -c

KMER.o: KMER.h KMER.cpp
	g++ -O2 -fopenmp KMER.h KMER.cpp -c

BaseStr.o: BaseStr.h BaseStr.cpp
	g++ -O2 -fopenmp BaseStr.h BaseStr.cpp -c

USet.o: USet.cpp USet.h
	g++ -O2 -fopenmp USet.cpp USet.h -c

ULLN.o: ULLN.h ULLN.cpp
	g++ -O2 -fopenmp ULLN.h ULLN.cpp -c
	
Utils.o: Utils.cpp Utils.h
	g++ -O2 -fopenmp Utils.cpp Utils.h -c

BWTs.o: BWTs.cpp BWTs.h
	g++ -O2 -fopenmp BWTs.cpp BWTs.h -c

BWTDtStr.o: BWTDtStr.cpp BWTDtStr.h
	g++ -O2 -fopenmp BWTDtStr.cpp BWTDtStr.h -c

GLOBAL.o: GLOBAL.cpp GLOBAL.h
	g++ -O2 GLOBAL.cpp GLOBAL.h -c

clean:
	rm -rf *.o *.gch MetaClusterTA

