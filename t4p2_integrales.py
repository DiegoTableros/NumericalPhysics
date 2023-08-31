#Tarea 4, problema 2. Calculo de integrales
#Diego Emilio Dominguez Tableros

import numpy as np
import matplotlib.pyplot as plt
from pylab import *

k=100000             #Numero de variables aleatorias
u=np.random.rand(k) #Variables aleatorias

#Estimacion de int_0^1 (1-x^2)^(3/2) dx
g=(1/k)*(1-u**2)**(3/2)     #Evaluacion de la funcion
res=sum(g)          #Promedio
print("Se estima que el valor de la primera integral es:\n",res)

#Estimacion de int_{-2}^2 exp(x+x^2) dx
h=(4/k)*exp((4*u-2)+(4*u-2)**2)
res=sum(h)
print("Se estima que el valor de la segunda integral es:\n",res)