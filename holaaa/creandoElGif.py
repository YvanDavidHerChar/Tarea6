import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

#Creamos una funcion que tome la gran lista y la divida en la mitad
def dividirLista(lista):
    mitad = len(lista)/2
    return lista[:mitad], lista[mitad:]

#Creamos 500 archivos .png a partir de los 500 archivos .dat 
for i in range(100):
    archivo = np.genfromtxt('laShit'+str(i*5)+'.dat')
    matrixU1, matrixV1 = dividirLista(archivo)
    matrixU = np.reshape(matrixU1, (41,41))
    matrixV = np.reshape(matrixV1, (41,41))

    nx = 41
    ny = 41

    x = np.linspace(0,2,nx)
    y = np.linspace(0,2,ny)

    fig = plt.figure(figsize=(11,7), dpi=100)
    ax = plt.gca(projection='3d')
    X,Y = np.meshgrid(x,y)
    wire1 = ax.plot_wireframe(X,Y,matrixU[:])
    wire2 = ax.plot_wireframe(X,Y,matrixV[:])
    plt.savefig('laVaina'+str(i)+'.png')
    plt.close(fig)

#Creamos un gif organizando los archivos para que no halla una discontinuidad en el gif.
os.system("convert -delay 10 $(ls -v *png) el2D.gif")
