#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float *matriz(int n_puntos)
{
  float *array;
  int i;
  int j;

  array = malloc(n_puntos*n_puntos*sizeof(float));

  for(i=0;i<n_puntos;i++)
  {
    for(j=0;j<n_puntos;j++)
     {
  	array[i+(j*n_puntos)] = 1.0;
    }
  }
  return array;
}

void print_matriz(float *matriz, int n_puntos)
{
  int i;
  int j;
  for(i=0 ; i<n_puntos ; i++)
  {
    for(j=0; j<n_puntos ; j++)
      {
	printf("%f\n", matriz[i + (j*n_puntos)]);
      }
  }
}

void copiar(float *origen, float *destino, int n_puntos)
{
  int i;
  int j;

  for(i=0; i < n_puntos ; i++)
      {
	for(j=0; j < n_puntos ; j++)
	  {
	    destino[i + (j*n_puntos)] = origen[i + (j*n_puntos)];
	  }
      }
}

void iteracion(float *u, float *v, float *un, float *vn, float dt, float dx, float dy, float nu, int n_puntos, int iteraciones)
{
  FILE *in;
  int t;
  int i;
  int j;
  int k;
  int l;
  char filename[100];

  for(t=0; t < iteraciones; t++)
    {
      copiar(u, un, n_puntos);
      copiar(v, vn, n_puntos);

      for(i=0; i<n_puntos; i++)
	{
	  for(j=0;j<n_puntos;j++)
	    {  
	      u[i + (j*n_puntos)] = un[i + (j*n_puntos)] - dt/dx*un[i + (j*n_puntos)]*(un[i + (j*n_puntos)]-un[(i-1) + (j*n_puntos)])-dt/dy*vn[i + (j*n_puntos)]*(un[i + (j*n_puntos)]-un[i + ((j-1)*n_puntos)])+nu*dt/(dx*dx)*(un[(i+1) + (j*n_puntos)]-2*un[i + (j*n_puntos)]+un[(i-1) + (j*n_puntos)])+nu*dt/(dy*dy)*(un[i + ((j+1)*n_puntos)]-2*un[i + (j*n_puntos)]+un[i + ((j+1)*n_puntos)]);
    
	      v[i + (j*n_puntos)] = vn[i + (j*n_puntos)] - dt/dx*vn[i + (j*n_puntos)]*(vn[i + (j*n_puntos)]-vn[(i-1) + (j*n_puntos)])-dt/dy*vn[i + (j*n_puntos)]*(vn[i + (j*n_puntos)]-vn[i + ((j-1)*n_puntos)])+nu*dt/(dx*dx)*(vn[(i+1) + (j*n_puntos)]-2*vn[i + (j*n_puntos)]+vn[(i-1) + (j*n_puntos)])+nu*dt/(dy*dy)*(vn[i + ((j+1)*n_puntos)]-2*vn[i + (j*n_puntos)]+vn[i + ((j+1)*n_puntos)]);
	    }

	  u[0 + (i*n_puntos)] = 1;
	  u[(n_puntos-1) + (i*n_puntos)] = 1;
	  u[i + (0*n_puntos)] = 1;
	  u[i + ((n_puntos-1)*n_puntos)] = 1;
    
	  v[0 + (i*n_puntos)] = 1;
	  v[(n_puntos-1) + (i*n_puntos)] = 1;
	  v[i + (0*n_puntos)] = 1;
	  v[i + ((n_puntos-1)*n_puntos)] = 1;
	}


  sprintf(filename, "laShit%d.dat", t);
  in = fopen(filename,"w"); 

  for(k=0 ; k<n_puntos ; k++)
  {
    for(l=0; l<n_puntos ; l++)
      {
	fprintf(in, "%f\n", u[k + (l*n_puntos)]);
      }
  }
  for(k=0 ; k<n_puntos ; k++)
  {
    for(l=0; l<n_puntos ; l++)
      {
	fprintf(in, "%f\n", v[k + (l*n_puntos)]);
      }
  }

  fclose(in);
    }
}

void estadoInicial(float *elU, float *elV, int nx)
{
  //Al instalar 41 puntos iniciales tenemos 11 puntos entre 0.5 y 1, pues estan incluidos los bordes.

  int puntos = 11;
  int i;
  int j;
  for(i=0; i < puntos; i++)
    {
      for(j=0; j < puntos; j++)
	{
	  elU[(10+i) + ((10+j)*nx)] = 2.0;
	  elV[(10+i) + ((10+j)*nx)] = 2.0;
	}
    }
}

int main()
{
  //inializando las constantes

  int nx, ny, iteraciones, c;
  float dx, dy, sigma, nu, dt;

  // inializando las matrices
  float *elU;
  float *elV;
  float *elOtroU;
  float *elOtroV;

  // Otorgando valores a las constantes

  nx = 41;
  ny = 41;
  iteraciones = 500;
  c = 1;
  dx = 2.0/(nx*1.0-1);
  dy = 2.0/(ny*1.0-1);
  sigma = .0009;
  nu = 0.01;
  dt = sigma*dx*dy/nu;

  elU = matriz(ny);
  elV = matriz(nx);
  elOtroU = matriz(nx);
  elOtroV = matriz(ny);

  estadoInicial(elU, elV, nx);
  iteracion(elU, elV, elOtroU, elOtroV, dt, dx, dy, nu,  nx, iteraciones);

  return 0;
}
