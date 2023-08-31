#Tarea 4, problema 3. Calculo de pi
#Diego Emilio Dominguez Tableros

import numpy as np

#APROXIMACION POR RAZON DE VOLUMENES
n=1000000       #Numeros aleatorios
a=-1            #Limite inferior
b=1             #Limite superior

dentro=0            #Numero de puntos dentro de esfera
for i in range(n):
    #Punto aleatorio en [0,1]^3
    p=np.random.rand(3)
    #Transformacion a [-1,1]^3
    p=a+(b-a)*p
    #Chequeo dentro esfera
    if(p[0]**2+p[1]**2+p[2]**2<=1):
        dentro+=1

pi_aprox=(6*dentro)/n
print("El calculo de pi por razon de volumenes por Montecarlo con ",n," puntos es de:")
print(pi_aprox)

#APROXIMACION POR INTEGRAL
k=1000000             #Numero de variables aleatorias
u=np.random.rand(k) #Variables aleatorias

#Estimacion de int_{-1}^1 (1-x^2)^(1/2) dx
g=(2/k)*np.sqrt(1-(2*u-1)**2)   #Evaluacion de la funcion
res=sum(g)                      #Promedio
#El resultado de la integral es pi/2
pi_aprox=2*res
print("El calculo de pi por aproximacion de integral por Montercalo con ",k," numeros aleatorios es de:")
print(pi_aprox)