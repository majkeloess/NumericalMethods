#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <cmath>
#include <fstream>

void fillMatrix(gsl_matrix *A)
{
  for (size_t i = 0; i < A->size1; ++i)
  {
    for (size_t j = 0; j < A->size2; ++j)
    {
      gsl_matrix_set(A, i, j, pow(2 + abs((int)i + (int)j), -abs((int)i + (int)j) / 2.0));
    }
  }
}

void powerMethod(gsl_matrix *A, int maxIter)
{
  size_t size = A->size1;
  gsl_vector *x = gsl_vector_alloc(size);
  gsl_vector_set_all(x, 1.0);
  double lambda = 0.0;

  std::ofstream file;
  file.open("data.csv");

  for (int i = 0; i < maxIter; ++i)
  {
    gsl_vector *Ax = gsl_vector_alloc(size);
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, Ax);
    double xAx;
    gsl_blas_ddot(x, Ax, &xAx);
    double xx;
    gsl_blas_ddot(x, x, &xx);
    lambda = xAx / xx;
    gsl_vector_memcpy(x, Ax);
    gsl_vector_scale(x, 1.0 / gsl_blas_dnrm2(x));

    file << i + 1 << '\t';
    for (size_t j = 0; j < size; ++j)
    {
      file << gsl_vector_get(x, j) << "\t";
    }
    file << std::endl;

    gsl_vector_free(Ax);
  }

  file.close();
  gsl_vector_free(x);
}

int main()
{
  int size = 7;
  gsl_matrix *A = gsl_matrix_alloc(size, size);
  fillMatrix(A);
  powerMethod(A, 120);
  gsl_matrix_free(A);
  return 0;
}
