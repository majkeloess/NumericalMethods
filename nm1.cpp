#include <fstream>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

double matrix_infinity_norm(gsl_matrix *A, int n)
{
  double max_sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    double row_sum = 0.0;
    for (int j = 0; j < n; j++)
    {
      row_sum += fabs(gsl_matrix_get(A, i, j));
    }
    if (row_sum > max_sum)
      max_sum = row_sum;
  }
  return max_sum;
}

void matrix_multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C, int n)
{
  // for (int i = 0; i < n; i++)
  // {
  //   for (int j = 0; j < n; j++)
  //   {
  //     double sum = 0;
  //     for (int k = 0; k < n; k++)
  //     {
  //       sum += gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j);
  //     }
  //     gsl_matrix_set(C, i, j, sum);
  //   }
  // }
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C);
}

int main()
{

  int n{4};
  int m{4};

  gsl_matrix *A{gsl_matrix_calloc(n, m)};
  gsl_matrix *A_copy{gsl_matrix_calloc(n, m)}; // Copy of A
  gsl_matrix *A_inv{gsl_matrix_calloc(n, m)};
  gsl_matrix *product{gsl_matrix_calloc(n, m)}; //  AA^-1
  gsl_vector *x{gsl_vector_calloc(m)};
  gsl_vector *b{gsl_vector_calloc(m)};
  int signum{1};
  gsl_permutation *p{gsl_permutation_calloc(m)};

  // 1, 2
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      double value = 1.0 / (i + j + 2);
      gsl_matrix_set(A, i, j, value);
      gsl_matrix_set(A_copy, i, j, value);
    }
  }

  gsl_linalg_LU_decomp(A, p, &signum);

  std::ofstream myfile;

  myfile.open("NumericalLab01.txt");

  myfile << "Diagonal A: ";

  double determinant{1};

  for (size_t i = 0; i < 4; i++)
  {
    determinant *= gsl_matrix_get(A, i, i);
    myfile << gsl_matrix_get(A, i, i) << " ";
  }
  myfile << std::endl
         << "Det(A): " << std::fixed << std::setprecision(13) << determinant << std::endl;

  // 3

  for (int i = 0; i < n; i++)
  {
    gsl_vector_set_zero(b);
    gsl_vector_set(b, i, 1.0);
    gsl_linalg_LU_solve(A, p, b, x);
    gsl_matrix_set_col(A_inv, i, x);
  }
  myfile << std::endl;
  myfile << "A^-1 :" << std::endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      myfile << std::setw(10) << std::fixed << std::setprecision(0) << gsl_matrix_get(A_inv, i, j) << " ";
    }
    myfile << std::endl;
  }

  // 4

  matrix_multiply(A_copy, A_inv, product, n); // AA^-1
  myfile << std::endl;
  myfile << "AA^-1 :" << std::endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      myfile << std::setw(20) << std::setprecision(17) << gsl_matrix_get(product, i, j) << " ";
    }
    myfile << std::endl;
  }

  // 5
  double norm_A = matrix_infinity_norm(A_copy, n);
  double norm_A_inv = matrix_infinity_norm(A_inv, n);

  double cond = norm_A * norm_A_inv;
  myfile << std::endl;
  myfile << "Condition number: " << std::setprecision(0) << cond << std::endl;

  myfile.close();

  gsl_matrix_free(A);
  gsl_matrix_free(A_copy);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_matrix_free(A_inv);
  gsl_matrix_free(product);
  gsl_permutation_free(p);

  return 0;
}
