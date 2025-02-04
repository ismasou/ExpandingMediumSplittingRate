CXFLAGS = -lfmt

all: Static Expanding ExpandingOp Opacity

Static:
	g++ ./src/Polynomial.cpp ./src/Static.cpp -o Static.out -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(CXFLAGS)

Opacity:
	g++ ./src/Polynomial.cpp ./src/Opacity.cpp -o Opacity.out -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(CXFLAGS)

Expanding:
	g++ ./src/Polynomial.cpp ./src/Expanding.cpp -o Expanding.out -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(CXFLAGS)

ExpandingOp:
	g++ ./src/Polynomial.cpp ./src/OE-Expanding.cpp -o OEExp.out -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native
