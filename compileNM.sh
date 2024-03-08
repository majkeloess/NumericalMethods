echo "Compiling and running"

param=$1

g++ "$param.cpp" -lgsl -lgslcblas -lm -o "$param"

./$param
