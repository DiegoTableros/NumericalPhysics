#Tarea 6, Fisica Numerica. Diego Emilio Dominguez Tableros
#Problema 1. Interpolacion de Lagrange

import numpy as np
import scipy as sc
from scipy.interpolate import lagrange
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

#Lectura de datos de archivo BretWigner.txt
datos=np.loadtxt("BretWigner.txt",skiprows=1,dtype=float)
E=datos[:,1]
f_E=datos[:,2]
sigma=datos[:,3]
n=len(E)

#Interpolacion de Lagrange con scipy.interpolate.lagrange
pol=lagrange(E,f_E)     #pol es una funcion evaluable

#Graficacion
n=int((E[n-1]-E[0])/1)          #Numero de puntos
x_aux=np.linspace(0,200,n)      #Arreglo para evaluar "pol"
f_aux=pol(x_aux)                #Evaluacion de "pol"
#Gráfica interpolacion
plt.plot(x_aux,f_aux,"g",label="Interp. Lagrange")

#plt.plot(E,f_E,"r.")            #Datos originales (puntos)
#Graficacion de datos originales con barras de error para f(E)
error_f_E=[sigma,sigma]
plt.errorbar(E,f_E,yerr=error_f_E,fmt='.',label="Datos")

plt.xlabel("Energía E [MeV]")
plt.ylabel("f(E) [MeV]")
plt.legend()
plt.grid(which="both")
plt.show()

h_media=max(f_aux)/2    #Altura media