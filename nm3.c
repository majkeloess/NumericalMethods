#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"
#include <stdio.h>
#include <math.h>

int save_vector_to_file(gsl_vector *b, const size_t n, char *name)
{
  FILE *f;
  f = fopen(name, "w");
  if (f == NULL)
  {
    perror("Nie udalo sie otworzyc pliku do zapisu");
    return 0;
  }

  double x=0.1;
  for (int i = 0; i < n; i++)
  {
    fprintf(f, "%f\t%5g\t\n",x, gsl_vector_get(b, i));
    x+=0.1;
  }
  fprintf(f, "\n");
  fflush(f);
  fclose(f);
  return 1;
}



gsl_vector *matrix_vector_multiply(gsl_matrix *A, gsl_vector *x, const int n)
{
  gsl_vector *c = gsl_vector_calloc(n);
  float sum = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      sum += gsl_matrix_get(A, i, j) * gsl_vector_get(x, j);
    }
    gsl_vector_set(c, i, sum);
    sum = 0;
  }
  return c;
}

double Lambda(float x)
{
  const float l1 = 0.4;
  const float d1 = 40;

  const float l2 = 0.2;
  const float d2 = 30;

  const float l3 = 0.1;
  const float d3 = 30;

  if (x <= d1)
  {
    return l1;
  }
  else if (x <= d1 + d2)
  {
    return l2;
  }
  else if (x <= d1 + d2 + d3)
  {
    return l3;
  }
  else
  {
    return 0;
  }
}

void fill_A(gsl_matrix *A, int n)
{
  float var = 100. / (n + 1);
  float dx05 = var / 2;
  float tmpLambda;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j)
      {
        tmpLambda = -Lambda(dx05 + var * i) - Lambda(dx05 + var * (i + 1));
        gsl_matrix_set(A, i, j, tmpLambda);
      }
      else if (i == j + 1)
      {
        tmpLambda = Lambda(dx05 + var * i);
        gsl_matrix_set(A, i, j, tmpLambda);
      }
      else if (i == j - 1)
      {
        tmpLambda = Lambda(dx05 + var * (j));
        gsl_matrix_set(A, i, j, tmpLambda);
      }
    }
  }
}

void fill_b(gsl_vector *b, int n)
{
  float var = (100. / (n + 1));
  float dx05 = var / 2;
  float tmpLambda;
  float Tl = 1000, TR = 100;

  for (int i = 0; i < n; i++)
  {
    if (i == 0)
    {
      tmpLambda = -Lambda(dx05 + var * i);
      // printf("%f\n",dx05 + var * i);
      // printf("%f\n",Tl*tmpLambda);
      // printf("%d",i);
      // printf("\n");
      gsl_vector_set(b, i, Tl * tmpLambda);
      // for(int j=0;j<n; j++)
      // printf("%.5f\t", gsl_vector_get(t, j));
    }
    else if (i == n - 1)
    {
      tmpLambda = -Lambda(dx05 + var * i);
      // printf("%f\n",TR*tmpLambda);
      gsl_vector_set(b, i, TR * tmpLambda);
    }
    else
      gsl_vector_set(b, i, 0);
  }
}

void result(gsl_matrix *A, gsl_vector *t, gsl_vector *b, const int n)
{
  gsl_vector *iloczyn = gsl_vector_calloc(n);
  gsl_vector *r_k = gsl_vector_calloc(n);
  gsl_vector *p = gsl_vector_calloc(n);
  float alfa = 0;
  float suma = 0;
  float mian = 0;
  float iloczyn_skal = 0;
  int licznik = 0;
    FILE *f;
  f = fopen("datanew.csv", "w");
  if (f == NULL)
  {
    perror("Nie udalo sie otworzyc pliku do zapisu");
    exit(1);
  }

  for (int x = 0; x < 60000; x++)
  {
    // iloczyn A * t_k
    for (int i = 0; i < n; i++)
    {
      suma = 0.0;
      for (int j = 0; j < n; j++)
      {
        suma += gsl_matrix_get(A, i, j) * gsl_vector_get(t, j);
      }
      gsl_vector_set(iloczyn, i, suma);
    }
    // r_k
    for (int i = 0; i < n; i++)
    {
      gsl_vector_set(r_k, i, gsl_vector_get(b, i) - gsl_vector_get(iloczyn, i));
    }


    // iloczyn r_kt * r_k // licznik w 2 równaniu
    iloczyn_skal = 0;
    for (int i = 0; i < n; i++)
      iloczyn_skal += gsl_vector_get(r_k, i) * gsl_vector_get(r_k, i);

    // sprawdzenie czy iloczyn nie wychozi poza zakres podany w zadaniu 10^-6
    if (sqrt(iloczyn_skal) < 0.000001)
    {
      return;
    }


    // iloczyn A*r_k
    for (int i = 0; i < n; ++i)
    {
      float s = 0;
      for (int j = 0; j < n; ++j)
      {
        s += gsl_matrix_get(A, i, j)*gsl_vector_get(r_k, j);
      }
      gsl_vector_set(p, i, s);
    }
    mian = 0;
    for (int i = 0; i < n; i++)
    {
      mian += gsl_vector_get(p, i) * gsl_vector_get(r_k, i);
    }

    alfa = (iloczyn_skal / mian);

    for (int i = 0; i < n; i++)
      gsl_vector_set(t, i, gsl_vector_get(t, i) + alfa * gsl_vector_get(r_k, i));

    float norma = 0;
    for (int i = 0; i < n; i++)
      norma += gsl_vector_get(t, i) * gsl_vector_get(t, i);
    fprintf(f,"%d\t%.2f\t%.2f\n", x, iloczyn_skal, sqrt(norma));
  }
  fclose(f);
}

int main()
{
  const int n = 9;


  gsl_matrix *A = gsl_matrix_calloc(n, n);
  gsl_vector *b = gsl_vector_calloc(n);
  gsl_vector *t = gsl_vector_calloc(n); // wyniki

  fill_A(A, n);
  fill_b(b, n);
  result(A, t, b, n);
  save_vector_to_file(t, n, "data.csv");

  return 0;
}