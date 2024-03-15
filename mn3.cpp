#include <iostream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

double Lambda(double x)
{
  if (x <= 40)
  {
    return 0.4;
  }
  else if (x <= 70)
  {
    return 0.2;
  }
  else
  {
    return 0.1;
  }
}

void wypelnij(gsl_matrix *A, int wymiar)
{
  double var = 100.0 / (wymiar + 1);
  double dx05 = var / 2;
  double tmpLambda;
  for (int i = 0; i < wymiar; i++)
  {
    for (int j = 0; j < wymiar; j++)
    {
      if (i == j)
      {
        tmpLambda = -Lambda(dx05 + var * (i + 1)) - Lambda(dx05 + var * i);
        gsl_matrix_set(A, i, j, tmpLambda);
      }
      else if (i == j + 1)
      {
        tmpLambda = Lambda(dx05 + var * i);
        gsl_matrix_set(A, i, j, tmpLambda);
      }
      else if (i == j - 1)
      {
        tmpLambda = Lambda(dx05 + var * j);
        gsl_matrix_set(A, i, j, tmpLambda);
      }
    }
  }
}

void wypelnij_wynik(gsl_vector *t, int wymiar)
{
  double var = 100.0 / (wymiar + 1);
  double dx05 = var / 2;
  double tmpLambda;
  double Tl = 1000, TR = 100;
  for (int i = 0; i < wymiar; i++)
  {
    if (i == 0)
    {
      tmpLambda = -Lambda(dx05 + var * i);
      gsl_vector_set(t, i, Tl * tmpLambda);
    }
    else if (i == wymiar - 1)
    {
      tmpLambda = -Lambda(dx05 + var * i);
      gsl_vector_set(t, i, TR * tmpLambda);
    }
    else
    {
      gsl_vector_set(t, i, 0);
    }
  }
}

void rozwiaz(gsl_matrix *a, gsl_vector *t, gsl_vector *wynik, int wymiar, std::ofstream &file)
{
  gsl_vector *iloczyn = gsl_vector_calloc(wymiar);
  gsl_vector *r = gsl_vector_calloc(wymiar);
  gsl_vector *p = gsl_vector_calloc(wymiar);
  double alfa = 0;
  double suma = 0.0;
  double mian = 0;
  double iloczyn_skal = 0;
  int licznik = 0;

  for (int x = 0; x < 60000; x++)
  {
    for (int i = 0; i < wymiar; i++)
    {
      suma = 0.0;
      for (int j = 0; j < wymiar; j++)
      {
        suma += gsl_matrix_get(a, i, j) * gsl_vector_get(t, j);
      }
      gsl_vector_set(iloczyn, i, suma);
    }

    for (int i = 0; i < wymiar; i++)
    {
      gsl_vector_set(r, i, gsl_vector_get(wynik, i) - gsl_vector_get(iloczyn, i));
    }

    iloczyn_skal = 0;
    for (int i = 0; i < wymiar; i++)
    {
      iloczyn_skal += gsl_vector_get(r, i) * gsl_vector_get(r, i);
    }

    if (sqrt(iloczyn_skal) < 0.000001)
    {
      return;
    }

    for (int i = 0; i < wymiar; ++i)
    {
      double s = 0;
      for (int j = 0; j < wymiar; ++j)
      {
        s += gsl_vector_get(r, j) * gsl_matrix_get(a, i, j);
      }
      gsl_vector_set(p, i, s);
    }

    mian = 0;
    for (int i = 0; i < wymiar; i++)
    {
      mian += gsl_vector_get(p, i) * gsl_vector_get(r, i);
    }

    alfa = (iloczyn_skal / mian);

    for (int i = 0; i < wymiar; i++)
    {
      gsl_vector_set(t, i, gsl_vector_get(t, i) + alfa * gsl_vector_get(r, i));
    }

    double norma = 0;
    for (int i = 0; i < wymiar; i++)
    {
      norma += gsl_vector_get(t, i) * gsl_vector_get(t, i);
    }

    file << std::fixed << x << " " << iloczyn_skal << " " << norma << std::endl;
  }
}

int main()
{

  int wymiar = 9; // lub 99
  gsl_matrix *A = gsl_matrix_calloc(wymiar, wymiar);
  gsl_vector *t = gsl_vector_calloc(wymiar);
  gsl_vector *wynik = gsl_vector_calloc(wymiar);

  wypelnij(A, wymiar);
  wypelnij_wynik(t, wymiar);

  std::ofstream myfile;
  myfile.open("data.csv");

  rozwiaz(A, t, wynik, wymiar, myfile);
  myfile << "Wyniki: \n";

  for (int i = 0; i < wymiar; i++)
  {
    myfile << std::setw(10) << std::setprecision(2) << gsl_vector_get(wynik, i) << std::endl;
  }

  gsl_matrix_free(A);
  gsl_vector_free(t);
  gsl_vector_free(wynik);

  return 0;
}
