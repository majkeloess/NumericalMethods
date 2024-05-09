#include <iostream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <random>
#include <chrono>

double fun(double x)
{
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);
  double delta = 5 * distribution(generator);
  return -0.25 * std::pow(x, 3) - 0.5 * std::pow(x, 2) + 5 * x + 5 + delta;
}

int main()
{
  for (int n : {11, 51, 101})
  {
    double x0 = -5, x1 = 5;
    gsl_matrix *X = gsl_matrix_calloc(n, 4);
    gsl_vector *Y = gsl_vector_calloc(n);
    int licznik = 0;
    for (double i = x0; i <= x1; i += (x1 - x0) / (n - 1))
    {
      for (int j = 0; j < 4; j++)
      {
        gsl_matrix_set(X, licznik, j, std::pow(i, j));
      }
      gsl_vector_set(Y, licznik, fun(i));
      licznik++;
    }

    gsl_matrix *XT = gsl_matrix_calloc(4, n);
    gsl_matrix_transpose_memcpy(XT, X);

    gsl_matrix *D = gsl_matrix_calloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, XT, X, 0.0, D);

    gsl_vector *r = gsl_vector_calloc(4);
    gsl_blas_dgemv(CblasNoTrans, 1.0, XT, Y, 0.0, r);

    gsl_vector *b = gsl_vector_calloc(4);
    gsl_permutation *p = gsl_permutation_calloc(4);
    int signum;
    gsl_linalg_LU_decomp(D, p, &signum);
    gsl_linalg_LU_solve(D, p, r, b);

    gsl_vector *Y_approx = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1.0, X, b, 0.0, Y_approx);

    std::ofstream file;
    file.open("data" + std::to_string(n) + ".csv");
    for (int i = 0; i < n; i++)
    {
      file << gsl_matrix_get(X, i, 1) << "\t" << gsl_vector_get(Y, i) << "\t" << gsl_vector_get(Y_approx, i) << std::endl;
    }
    file.close();

    gsl_matrix_free(X);
    gsl_matrix_free(XT);
    gsl_matrix_free(D);
    gsl_vector_free(Y);
    gsl_vector_free(r);
    gsl_vector_free(b);
    gsl_vector_free(Y_approx);
    gsl_permutation_free(p);
  }
  return 0;
}
