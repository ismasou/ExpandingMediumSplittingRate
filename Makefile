FMT_FLAGS = -I./external/fmt/include -L./external/fmt/build -lfmt
GSL_INC   = -I$(HOME)/software/noctua2/gsl/include
GSL_LIB   = -L$(HOME)/software/noctua2/gsl/lib

all: Static Static2D Expanding ExpandingOp Opacity

Static:
	g++ $(GSL_INC) ./src/Polynomial.cpp ./src/Static.cpp -o Static.out $(GSL_LIB) -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(FMT_FLAGS)

Static2D:
	g++ $(GSL_INC) ./src/Polynomial.cpp ./src/Static2DIntegral.cpp -o Static2D.out $(GSL_LIB) -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(FMT_FLAGS)

TestTheta:
	g++ ./src/TestThetaIntegral.cpp -o TestTheta.out -lm -O3 $(FMT_FLAGS)

Opacity:
	g++ $(GSL_INC) ./src/Polynomial.cpp ./src/Opacity.cpp -o Opacity.out $(GSL_LIB) -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(FMT_FLAGS)

Expanding:
	g++ $(GSL_INC) ./src/Polynomial.cpp ./src/Expanding.cpp -o Expanding.out $(GSL_LIB) -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(FMT_FLAGS)

ExpandingOp:
	g++ $(GSL_INC) ./src/Polynomial.cpp ./src/OE-Expanding.cpp -o OEExp.out $(GSL_LIB) -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(FMT_FLAGS)

CheckInterp:
	g++ $(GSL_INC) ./src/check_interp.cpp -o check_interp.out $(GSL_LIB) -lm -lgsl -lgslcblas -g -O3 -fopenmp -march=native $(FMT_FLAGS)
