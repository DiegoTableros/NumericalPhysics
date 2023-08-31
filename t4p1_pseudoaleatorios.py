#Tarea 4, problema 1. Generacion de num pseudoaleatorios
#Diego Emilio Dominguez Tableros

import numpy as np
import matplotlib.pyplot as plt

n=10                #Cantidad de aleatorios
ale=np.zeros(n)     #Arreglo para numeros
ale[0]=semilla=12   #Valor inicial (seed)
a=57                #Parametro a
c=1                 #Parametro c
m=256               #Parametro M

#Generacion de algunos numeros pseudoaleatorios
for i in range(n):
    ale[i]=(a*ale[i-1]+c)%m
print("Algunos numeros pseudoaleatorios generados por el método de congruencias lineales con parametros")
print("a=",a,", c=",c,", M=",m,", seed=",semilla)
print(ale)
print("\n")

#Determinacion de periodo
ale=semilla=10
new_ale=(a*ale+c)%m
i=1
while(new_ale!=semilla):
    ale=new_ale
    new_ale=(a*ale+c)%m
    i+=1
print("El periodo del método de congruencias lineales con parametros a=",a,", c=",c,", M=",m,", seed=",semilla,"es de:\n",i,"\n")

#Generación de lista de aleatorios
ale=np.zeros(i)
ale[0]=semilla=10
for i in range(len(ale)):
    ale[i]=(a*ale[i-1]+c)%m

#Graficacion (x_{2i-1},x_{2i}), i=1,2,...
aux=int(len(ale)/2)
plt.figure(1)
for i in range(aux):
    plt.plot(ale[2*i-1],ale[2*i],'.')
plt.xlabel("x_{2i-1}")
plt.ylabel("y_{2i}")

#Graficacion (x_i,i)
aux=int(len(ale))
plt.figure(2)
for i in range(aux):
    plt.plot(i,ale[i],'.')
plt.xlabel("i")
plt.ylabel("x_i")