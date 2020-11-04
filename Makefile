all: test test2

FLAGS=-std=c++17 -Wall -g -O0

DEFS=-DNO_AFFIN_EXPR -DEXPLICIT_TEMPLATES

LIBGEOM=/home/salvi/project/libgeom

SKETCHES=/home/salvi/project/sketches/engine

INCLUDES=-I/usr/include/eigen3 \
	-I$(SKETCHES)/src/Affin \
	-I$(SKETCHES)/src/BaseFunctions \
	-I$(SKETCHES)/src/Container \
	-I$(SKETCHES)/src/Curves \
	-I$(SKETCHES)/src/Globals \
	-I$(SKETCHES)/src/Mesh \
	-I$(SKETCHES)/src/Numeric/Misc \
	-I$(SKETCHES)/src/Surfaces

LIBS=\
	-L$(SKETCHES)/debug/Affin -lAffin \
	-L$(SKETCHES)/debug/BaseFunctions -lBaseFunctions \
	-L$(SKETCHES)/debug/Container -lContainer \
	-L$(SKETCHES)/debug/Curves -lCurves \
	-L$(SKETCHES)/debug/Globals -lGlobals \
	-L$(SKETCHES)/debug/Numeric/LinAlgTypes -lLinAlgTypes \
	-L$(SKETCHES)/debug/Numeric/Misc -lMisc \
	-L$(SKETCHES)/debug/Surfaces -lSurfaces

test: test.cc
	$(CXX) -o $@ $< $(FLAGS) $(DEFS) $(INCLUDES) $(LIBS) $(LIBS)

test2: test2.cc
	$(CXX) -o $@ $< $(FLAGS) -I/usr/include/eigen3 -I$(LIBGEOM) -L$(LIBGEOM)/debug -lgeom
