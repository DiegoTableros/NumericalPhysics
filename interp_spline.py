#Tarea 6, Fisica Numerica. Diego Emilio Dominguez Tableros
#Problema 2. Interpolacion por Splines cubicos

import numpy as np
import scipy as sc
from scipy.interpolate import CubicSpline
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

#Lectura de datos de archivo BretWigner.txt
datos=np.loadtxt("BretWigner.txt",skiprows=1,dtype=float)
E=datos[:,1]
f_E=datos[:,2]
sigma=datos[:,3]
n=len(E)

#Interpolacion de Lagrange con scipy.interpolate.lagrange
spline=CubicSpline(E,f_E)     #spline es una funcion evaluable

#Graficacion
n=int((E[n-1]-E[0])/1)          #Numero de puntos
x_aux=np.linspace(0,200,n)      #Arreglo para evaluar "spline"
f_aux=spline(x_aux)             #Evaluacion de "spline"
#Gráfica interpolacion
plt.plot(x_aux,f_aux,"y",label="Interp. Splines cúbicos")

#plt.plot(E,f_E,"r.")            #Datos originales (puntos)
#Graficacion de datos originales con barras de error para f(E)
error_f_E=[sigma,sigma]
plt.errorbar(E,f_E,yerr=error_f_E,fmt='.',label="Datos")

plt.xlabel("Energía E [MeV]")
plt.ylabel("f(E) [MeV]")
plt.legend()
plt.grid(which="both")
plt.show()