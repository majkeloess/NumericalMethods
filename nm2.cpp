#include <fstream>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void gsl_matrix_scale_row(gsl_matrix *A, int i, double factor)
{
  for (int j = 0; j < A->size2; j++)
  {
    double Aij = gsl_matrix_get(A, i, j);
    gsl_matrix_set(A, i, j, Aij * factor);
  }
}

void gsl_matrix_sub_row(gsl_matrix *A, int i, const gsl_matrix *B, int j, double factor)
{
  for (int k = 0; k < A->size2; k++)
  {
    double Aik = gsl_matrix_get(A, i, k);
    double Bjk = gsl_matrix_get(B, j, k);
    gsl_matrix_set(A, i, k, Aik - factor * Bjk);
  }
}

void gauss_jordan(gsl_matrix *A, gsl_vector *b)
{
  int n = A->size1;

  for (int i = 0; i < n; i++)
  {

    double Aii = gsl_matrix_get(A, i, i);
    gsl_matrix_scale_row(A, i, 1.0 / Aii);
    gsl_vector_set(b, i, gsl_vector_get(b, i) / Aii);

    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        double factor = gsl_matrix_get(A, j, i);
        gsl_matrix_sub_row(A, j, A, i, factor);
        gsl_vector_set(b, j, gsl_vector_get(b, j) - factor * gsl_vector_get(b, i));
      }
    }
  }
}

int main()
{

  std::ofstream myfile;
  myfile.open("NumericalLab02.txt");

  double A_data[5][5] = {
      {0, 1, 6, 9, 10},
      {2, 1, 6, 9, 10},
      {1, 6, 6, 8, 6},
      {5, 9, 10, 7, 10},
      {3, 4, 9, 7, 9}};
  double b_data[5] = {10, 2, 9, 9, 3};

  for (double q = 0; q <= 3; q += 0.01)
  {
    if (q != 1.00)
    {

      gsl_matrix *A = gsl_matrix_alloc(5, 5);
      gsl_vector *b = gsl_vector_alloc(5);
      for (int i = 0; i < 5; i++)
      {
        for (int j = 0; j < 5; j++)
        {
          gsl_matrix_set(A, i, j, A_data[i][j]);
        }
        gsl_vector_set(b, i, b_data[i]);
      }

      gsl_matrix_set(A, 0, 0, q);

      gauss_jordan(A, b);

      gsl_vector *c = gsl_vector_alloc(5);
      for (int i = 0; i < 5; i++)
      {
        double sum = 0;
        for (int j = 0; j < 5; j++)
        {
          sum += gsl_matrix_get(A, i, j) * gsl_vector_get(b, j);
        }
        gsl_vector_set(c, i, sum);
      }

      myfile << "q = " << std::fixed << std::setw(5) << q << " , x = [";
      for (int i = 0; i < 5; i++)
      {
        myfile << gsl_vector_get(b, i) << std::setw(5);
        if (i < 4)
          myfile << ", ";
      }
      myfile << "] ";

      myfile << " , c = [";
      for (int i = 0; i < 5; i++)
      {
        myfile << gsl_vector_get(c, i) << std::setw(5);
        if (i < 4)
          myfile << ", ";
      }
      myfile << "]\n";

      gsl_matrix_free(A);
      gsl_vector_free(b);
      gsl_vector_free(c);
    }
  }

  myfile.close();

  return 0;
}