MetaClusterTA:  MetaClusterTA.cpp Methods.o MetaCluster.o TomAlgorithm.o BaseStr.o ULLN.o KMER.o ShortKMER.o USet.o Reader.o Utils.o GLOBAL.o BWTs.o BWTDtStr.o Structs.o
	g++ -fopenmp MetaClusterTA.cpp Methods.o MetaCluster.o TomAlgorithm.o BaseStr.o ULLN.o KMER.o ShortKMER.o USet.o Reader.o Utils.o GLOBAL.o BWTs.o BWTDtStr.o Structs.o -o MetaClusterTA

Methods.o: Methods.h Methods.cpp
	g++ -fopenmp Methods.h Methods.cpp -c

MetaCluster.o: MetaCluster.cpp MetaCluster.h
	g++ -fopenmp MetaCluster.cpp MetaCluster.h -c

Reader.o: Reader.h Reader.cpp
	g++ -fopenmp Reader.h Reader.cpp -c

TomAlgorithm.o: TomAlgorithm.cpp TomAlgorithm.h
	g++ -fopenmp TomAlgorithm.cpp  TomAlgorithm.h -c

ShortKMER.o: ShortKMER.h ShortKMER.cpp
	g++ -fopenmp ShortKMER.h ShortKMER.cpp -c

KMER.o: KMER.h KMER.cpp
	g++ -fopenmp KMER.h KMER.cpp -c

BaseStr.o: BaseStr.h BaseStr.cpp
	g++ -fopenmp BaseStr.h BaseStr.cpp -c

USet.o: USet.cpp USet.h
	g++ -fopenmp USet.cpp USet.h -c

ULLN.o: ULLN.h ULLN.cpp
	g++ -fopenmp ULLN.h ULLN.cpp -c
	
Utils.o: Utils.cpp Utils.h
	g++ -fopenmp Utils.cpp Utils.h -c

BWTs.o: BWTs.cpp BWTs.h
	g++ -fopenmp BWTs.cpp BWTs.h -c

BWTDtStr.o: BWTDtStr.cpp BWTDtStr.h
	g++ -fopenmp BWTDtStr.cpp BWTDtStr.h -c

Structs.o: Structs.cpp Structs.h
	g++ -fopenmp Structs.cpp Structs.h -c

GLOBAL.o: GLOBAL.cpp GLOBAL.h
	g++ GLOBAL.cpp GLOBAL.h -c

clean:
	rm -rf *.o *.gch MetaClusterTA

