#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define pi  3.141592653589793
#define c   300000000.0
#define Rt  6378100.0
#define m   (1.67262177*pow(10.,-27.))
#define e   (1.602176287*pow(10., -19))
#define Bo  (3.*pow(10.,-5.))

float* crearVector(float v1, float v2, float v3)
{
  float* vector;

  vector = malloc(3*sizeof(float));

  vector[0] = v1;
  vector[1] = v2;
  vector[2] = v3;

  return vector;
}


float  magnitud( float* vector)
{
  int i;
  float suma = 0;
  float resultado = 0;

  suma =  pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2);

  resultado = sqrt(suma);

  return resultado; 
}


void copy(float *origen, float *destino, int n_puntos){
  int i;
  for(i=0;i<n_puntos;i++){
    destino[i] = origen[i];
  }
}

float* calcularAceleracion(float* pos , float* vel,float gamma)
{
  float* acel;
  float constante;
  float r;
  int i;

  acel = crearVector(0,0,0);
  r = magnitud(pos);
  constante = -((Bo*pow(Rt,3)*e)/(pow(r,5)*gamma*m));
  acel[0] = constante*(vel[1]*(2*pos[2]*pos[2] - pos[0]*pos[0] - pos[1]*pos[1]) - 3*vel[2]*pos[0]*pos[2]);
  acel[1] = -1*constante*(vel[0]*(2*pos[2]*pos[2] - pos[0]*pos[0] - pos[1]*pos[1]) - 3*vel[2]*pos[0]*pos[2]);
  acel[2] = constante*(vel[0]*3*pos[2]*pos[1] - vel[1]*3*pos[0]*pos[2]);

  return acel;
}

void RungeKutta(float* posI , float* velI, int iteraciones, float paso, float energia, float angulo, float gamma)
{
  //Creamos e inicializamos todas las variables
  FILE *in;
  int i;
  float* elVectorX;
  float* elVectorV;
  float* elVectorFinalX;
  float* elVectorFinalV;
  float* aceleracion1;
  float* aceleracion2;
  float* aceleracion3;
  float* aceleracion4;
  float* vel2;
  float* vel3;
  float* vel4;
  float* pos2;
  float* pos3;
  float* pos4;
  float k1x, k1y, k1z , k1vx, k1vy, k1vz;
  char filename[100];

  elVectorX = posI;
  elVectorV = velI;
  elVectorFinalX = crearVector(0,0,0);
  elVectorFinalV = crearVector(0,0,0);
  aceleracion1 = crearVector(0,0,0);
  aceleracion2 = crearVector(0,0,0);
  aceleracion3 = crearVector(0,0,0);
  aceleracion4 = crearVector(0,0,0);
  vel2 = crearVector(0,0,0);
  vel3 = crearVector(0,0,0);
  vel4 = crearVector(0,0,0);
  pos2 = crearVector(0,0,0);
  pos3 = crearVector(0,0,0);
  pos4 = crearVector(0,0,0);

  //Creamos el archivo en el que vamos a escribir la informacion
  sprintf(filename, "trayectoria_%d_%d.dat",(int) energia,(int) angulo*180/pi);
  in = fopen(filename,"w");   
  fprintf(in, "%d %f %f %f\n", 0, elVectorX[0], elVectorX[1], elVectorX[2]); 
  fclose(in);

  for(i=0 ; i<iteraciones ; i++)
    {
      //Comienzo
      k1x = elVectorX[0];
      k1y = elVectorX[1];
      k1z = elVectorX[2];
      k1vx = elVectorV[0];
      k1vy = elVectorV[1];
      k1vz = elVectorV[2];
      aceleracion1 = calcularAceleracion(elVectorX, elVectorV, gamma);
  
    //Segundo Paso
      pos2[0] = elVectorX[0] + 0.5*paso*k1vx;
      pos2[1] = elVectorX[1] + 0.5*paso*k1vy;
      pos2[2] = elVectorX[2] + 0.5*paso*k1vz;
      vel2[0] = elVectorV[0] + 0.5*paso*aceleracion1[0];
      vel2[1] = elVectorV[1] + 0.5*paso*aceleracion1[1];
      vel2[2] = elVectorV[2] + 0.5*paso*aceleracion1[2];
      aceleracion2 = calcularAceleracion(pos2, vel2, gamma);

      //Tercer Paso
      pos3[0] = elVectorX[0] + 0.5*paso*vel2[0];
      pos3[1] = elVectorX[1] + 0.5*paso*vel2[1];
      pos3[2] = elVectorX[2] + 0.5*paso*vel2[2];
      vel3[0] = elVectorV[0] + 0.5*paso*aceleracion2[0];
      vel3[1] = elVectorV[1] + 0.5*paso*aceleracion2[1];
      vel3[2] = elVectorV[2] + 0.5*paso*aceleracion2[2];
      aceleracion3 = calcularAceleracion(pos3,vel3, gamma);

      //Cuarto Paso
      pos4[0] = elVectorX[0] + paso*vel3[0];
      pos4[1] = elVectorX[1] + paso*vel3[1];
      pos4[2] = elVectorX[2] + paso*vel3[2];
      vel4[0] = elVectorV[0] + paso*aceleracion3[0];
      vel4[1] = elVectorV[1] + paso*aceleracion3[1];
      vel4[2] = elVectorV[2] + paso*aceleracion3[2];
      aceleracion4 = calcularAceleracion(pos4,vel4, gamma);

      //Hacemos el Promedio
      elVectorFinalX[0] = elVectorX[0] + paso*(k1x + 2*pos2[0] + 2*pos3[0] + pos4[0])/6;
      elVectorFinalX[1] = elVectorX[1] + paso*(k1y + 2*pos2[1] + 2*pos3[1] + pos4[1])/6;
      elVectorFinalX[2] = elVectorX[2] + paso*(k1z + 2*pos2[2] + 2*pos3[2] + pos4[2])/6;
      elVectorFinalV[0] = elVectorV[0] + paso*(k1vx + 2*vel2[0] + 2*vel3[0] + vel4[0])/6;
      elVectorFinalV[1] = elVectorV[1] + paso*(k1vy + 2*vel2[1] + 2*vel3[1] + vel4[1])/6;
      elVectorFinalV[2] = elVectorV[2] + paso*(k1vz + 2*vel2[2] + 2*vel3[2] + vel4[2])/6;
     
      //Escribimos en el archivo final
      in = fopen(filename,"a");   
      fprintf(in, "%d %f %f %f\n", (i+1), elVectorFinalX[0], elVectorFinalX[1], elVectorFinalX[2]); 
      fclose(in);
      //Copiamos para volver a iterar
      copy(elVectorFinalX , elVectorX,3);
      copy(elVectorFinalV, elVectorV,3);
    }

}

int main(int argc, char **argv)
{
  float* posI;
  float* velI;
  float paso;
  int numeroDePuntos;
  float angulo;
  float energia;
  float nuevaEnergia;
  float v;
  float gamma;
 
 
  energia = atof(argv[1]);
  angulo = atof(argv[2]); 
  angulo = angulo*pi/180;
  nuevaEnergia = energia*e*pow(10.0,6);
  gamma = nuevaEnergia/(m*c*c)+1 ; 
  paso = 2*pi*gamma*m/(50*e*Bo);
  numeroDePuntos = 300000;
  v = sqrt(c*c*(1-(((m*c*c)/(nuevaEnergia+m*c*c))*((m*c*c)/(nuevaEnergia+m*c*c)))));
  posI = crearVector(2*Rt,0,0);
  velI = crearVector(0,v*sin(angulo), v*cos(angulo));

  RungeKutta(posI , velI, numeroDePuntos, paso, energia, angulo, gamma);
  return 0;
}
