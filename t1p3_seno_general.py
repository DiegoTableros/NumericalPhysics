#Tarea 1, problema 3. Cálculo de sen(x) general
#Diego Emilio Dominguez Tableros

from math import *
import time

#Valor para el cálculo de sen(x)
x=6*pi/5
print("Valor de x original: ",x,"\n");
#Escalamiento al intervalo x<2*pi
if x>(2*pi):
    v=int(x/(2*pi))                         #Cuántas veces cabe 2*pi en x
    x=x-(v*2*pi)                            #Quitamos "v" veces 2*pi del x
    print("Valor de x escalado: ",x,"\n");

#Relación de recurrencia del seno
def a_n(n,a_ant):
    #EL siguiente término "a_sig" es múltiplo del
    #término anterior "a_ant" y utiliza el número de
    #término "n" y el valor "x"
    a_sig=(-1)*((x**2)/((2*n-2)*(2*n-1)))*a_ant
    return a_sig

#Suma parcial de términos de la recurrencia del seno
def suma(N):
    a_act=x;            #El primer término de la recurrencia es el mismo "x"
    res=a_act;          #El primer sumando es el a_1
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
tol=5.0e-17
#Encabezado de tabla
print("%-6s%-24s%-24s" % ("N","Suma(N)","Error relativo"))
error_rel=1             #Valor inicial del error relativo del 100%
error_abs=1             #Valor inicial del error absoluto de 1
N=1                     #N inicial
#Mientras el error absoluto calculado sea mayor
#que la tolerancia, mejoramos nuestra aproximación
while error_abs>tol:
    #Calcula sen(x) con la función "suma(N)" para el N actual
    suma_N=suma(N)
    #Calculo del error relativo y absoluto con la función de la biblioteca math de Python
    error_abs=abs(suma_N-sin(x))
    error_rel=abs(error_abs/sin(x))
    #Impresión de resultados para el N actual
    print('%-6i%-24.16g%-24.16g' % (N,suma_N,error_rel))
    #time.sleep(0.5)
    #Incremento de N en 1 para la sig iteración
    N+=1