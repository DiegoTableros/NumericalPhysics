#Tarea 1, problema 3. Cálculo de sen(x)
#Diego Emilio Dominguez Tableros

from math import *

#Valor para el cálculo de sen(x)
x=pi/2
print("Valor de x: ",x,"\n");

#Relación de recurrencia del seno
def a_n(n,a_ant):
    #EL siguiente término "a_sig" es múltiplo del
    #término anterior "a_ant" y utiliza el número de
    #término "n" y el valor "x"
    a_sig=(-1)*((x**2)/((2*n-2)*(2*n-1)))*a_ant
    #Se retorna el término a_n
    return a_sig

#Suma parcial de términos de la recurrencia del seno
def suma(N):
    #El primer término de la recurrencia es el mismo "x"
    a_act=x;
    #El primer sumando es el a_1
    res=a_act;
    #Se calcula la suma a partir del a_2 hasta a_N
    for i in range(2,N+1):
        #Se actualiza al término actual a_i usando
        #la función "a_n()" con el término anterior
        a_act=a_n(i,a_act)
        #Se suma al resultado
        res+=a_act
    #Retorno de la suma desde a_1 hasta a_N
    return res

#Definición de la tolerancia
tol=1.0e-8
#Encabezado de tabla
print('%-3s%-20s%-20s' % ("N","Suma(N)","Error relativo"))
#Valor inicial del error relativo del 100%
error_rel=1
#N inicial
N=1
#Mientras el error relativo calculado sea mayor
#que la tolerancia, mejoramos nuestra aproximación
while error_rel>tol:
    #Calcula sen(x) con la función "suma(N)"
    #para el N actual
    suma_N=suma(N)
    #Calculo del error relativo con la función
    #de la biblioteca math de Python
    error_rel=abs((suma_N-sin(x))/sin(x))
    #Impresión de resultados para el N actual
    print('%-3i%-20.16f%-20.16g' % (N,suma_N,error_rel))
    #Incremento de N en 1 para la sig iteración
    N+=1