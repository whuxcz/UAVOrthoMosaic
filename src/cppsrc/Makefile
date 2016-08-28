######################################################################

KEYMATCHFULL=KeyMatchFull
GENDEM=genDem
MERGEPLY=mergePly
MODKEY=modKey
RANSAC=RANSAC.o
RANSACModel=RANSACModel.o
KEYS2A_FLANN=keys2a_flann.o
EXIF=exif.o
MULTIORTHO=Multiple
MULTIMOSAIC=Multimosaic

#######################################################################
CXX=g++

#Your INCLUDE path(e.g., -I/usr/include)
INCLUDE=-I/usr/include -I/usr/include/eigen3

#Your LIBRARY PATH(e.g.,-L /usr/lib/i386-linux-gnu -lboost_system)
LIB=-L/usr/lib/i386-linux-gnu -L/usr/local/lib

LDFLAG= -lz -lboost_system -lflann -ltiff -lgeotiff -ljpeg

#Complie Options
STD = -std=gnu++0x


########################################################################
all: GENDEM KEYMATCHFULL MERGEPLY MODKEY MULTIORTHO MULTIMOSAIC

GENDEM: RANSAC RANSACModel EXIF genDem.cpp
	$(CXX) -o $(GENDEM) genDem.cpp $(RANSAC) $(RANSACModel) $(EXIF) $(INCLUDE) $(LIB) $(LDFLAG) $(STD)
RANSAC: RANSACModel RANSAC.cpp RANSAC.h
	$(CXX) -c  RANSAC.cpp $(INCLUDE) $(LIB) $(LDFLAG) $(STD)
RANSACModel: RANSACModel.cpp RANSACModel.h
	$(CXX) -c  RANSACModel.cpp $(INCLUDE) $(LIB) $(LDFLAG) $(STD)
EXIF: exif.cpp exif.h
	$(CXX) -c exif.cpp $(INCLUDE) $(LIB) $(LDFLAG) $(STD)


KEYMATCHFULL: KEYS2A_FLANN KeyMatchFull_flann.cpp
	$(CXX) -o $(KEYMATCHFULL) KeyMatchFull_flann.cpp $(KEYS2A_FLANN) $(INCLUDE) $(LIB) $(LDFLAG) $(STD)
KEYS2A_FLANN: keys2a_flann.h keys2a_flann.cpp
	$(CXX) -c keys2a_flann.cpp $(INCLUDE_FLANN) $(LIB) $(LDFLAG) $(STD)

MERGEPLY: mergePly.cpp
	$(CXX) -o $(MERGEPLY) mergePly.cpp

MODKEY: modKey.cpp
	$(CXX) -o $(MODKEY) modKey.cpp

MULTIMOSAIC: h_multimosaic.cpp matrix.h imgutil.hpp
	$(CXX) -o $(MULTIMOSAIC) h_multimosaic.cpp $(INCLUDE) $(LIB) $(LDFLAG)

MULTIORTHO: h_multiple.cpp matrix.h imgutil.hpp
	$(CXX) -o $(MULTIORTHO) h_multiple.cpp $(INCLUDE) $(LIB) $(LDFLAG)

clean:
	rm -f *.o $(GENDEM) $(KEYMATCHFULL) $(MERGEPLY) $(MULTIORTHO) $(MULTIMOSAIC) $(MODKEY)
