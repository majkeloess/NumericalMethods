echo "Compiling and running"

param=$1

gcc "$param.c" -lgsl -lgslcblas -lm -o "$param"

./$param
