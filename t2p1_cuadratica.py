#Tarea 2, problema 1. Formula cuadratica
#Diego Emilio Dominguez Tableros
from numpy import sqrt

#Coeficientes a,b,c de la ec. ax²+bx+c=0
a=1.0
b=1.0
#Definicion de n maximo para c=10^(-n)
n_max=18;

for n in range(1,n_max+1):
    #Definicion del coeficiente c
    c=10**(-n)
    #Calculo de soluciones
    x_1=(-b+sqrt(b**2-4*a*c))/(2*a)
    x_2=(-b-sqrt(b**2-4*a*c))/(2*a)
    x_1p=(-2*c)/(b+sqrt(b**2-4*a*c))
    x_2p=(-2*c)/(b-sqrt(b**2-4*a*c))
    #Impresion de resultados
    print("Soluciones de ax²+bx+c=0 con a=",a,", b=",b,", c=",c)
    print('%-5s%-24.16g%-5s%-24.16g%-7s%-24.16g' % ("x_1=",x_1,"x'_1=",x_1p,"ErrAbs=",abs(x_1-x_1p)))
    print('%-5s%-24.16g%-5s%-24.16g%-7s%-24.16g' % ("x_2=",x_2,"x'_2=",x_2p,"ErrAbs=",abs(x_2-x_2p)))
    print("\n")