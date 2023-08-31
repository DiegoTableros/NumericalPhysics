#Tarea 2, problema 2. Funciones de Bessel esfericas
#Diego Emilio Dominguez Tableros

from numpy import sin
from numpy import cos
from numpy import zeros

x_val=[0.1,1,10]                #Valores de x
l_max=25                        #Maximo valor de l a calcular
sol=zeros((l_max+1,2),float)    #Matriz de soluciones
UP=0; DOWN=1;                   #Indicadores auxiliares
#Tabla de valores de j_3, j_5 y j_8 para x={0.1,1,10}
vals=zeros((4,4),float)

def hacia_arriba(x):
    #Primeros valores de j_l(x) ya son conocidos
    sol[0,UP]=sin(x)/x                      #j_0(x)
    sol[1,UP]=(sin(x)/(x**2))-(cos(x)/x)    #j_1(x)
    #Cálculo de j_l(x) hasta l_max
    for l in range(1,l_max):
        sol[l+1,UP]=((2*l+1)/x)*sol[l,UP]-sol[l-1,UP]

def hacia_abajo(x):
    #Definicion del L para comenzar arbitrariamente
    L=l_max
    #Valores arbitrarios
    j_lM1=0         #j_(L+1)
    j_l=1           #j_L
    #Calculo de j_l(x) hacia abajo
    for l in range(L,0,-1):
        #Recurrencia DOWN
        j_lm1=((2*l+1)/x)*j_l-j_lM1;
        #Guardado de valores que SI interesan (0<=l<=l_max)
        if (l-1)<=l_max:
            sol[l-1,DOWN]=j_lm1
        #Intercambio de valores para sig iteracion
        j_lM1=j_l
        j_l=j_lm1
    #En este punto ya se tiene j_0 computado
    j_0=sin(x)/x            #Valor real de j_0
    j_0C=sol[0,DOWN]        #Valor de j_0 Computado
    factor=j_0/j_0C         #Factor de normalizacion
    #Ya se pueden normalizar los valores
    for l in range(0,l_max+1):
        #j_l Normalizado = j_l Computado * Factor
        sol[l,DOWN]=sol[l,DOWN]*factor

vals[1,1]=9.518519719e-6
vals[1,2]=9.616310231e-10
vals[1,3]=2.901200102e-16
vals[2,1]=9.006581118e-3
vals[2,2]=9.256115862e-5
vals[2,3]=2.826498802e-8
vals[3,1]=-3.949584498e-2
vals[3,2]=-5.553451162e-2
vals[3,3]=1.255780236e-1

i=1            #Contador de valor de x actual
for x in x_val:
    #Cálculo de j_l(x) para el x actual
    hacia_arriba(x)
    hacia_abajo(x)
    #Impresion de resultados
    print("Valores de la función de Bessel esférica j_l(x) con x=",x)
    print('%-6s%-24s%-24s%-24s' % ("l","j_l UP","j_l DOWN","|UP-DOWN|/(|UP|+|DOWN|)"))
    for l in range(0,l_max+1):
        j_lUP=sol[l,UP]
        j_lDOWN=sol[l,DOWN]
        dif=abs(j_lUP-j_lDOWN)/(abs(j_lUP)+abs(j_lDOWN))
        print('%-6i%-24.16g%-24.16g%-24.16g' % (l,j_lUP,j_lDOWN,dif))
    print("\n")
    
    j=1         #Contador auxiliar
    print("Comparación con valores dados")
    print('%-8s%-30s%-30s' % ("j_l","Error relativo con j_l UP","Error relativo con j_l DOWN"))
    for n in [3,5,8]:
        errUP=abs(sol[n,UP]-vals[i,j])/abs(vals[i,j])
        errDOWN=abs(sol[n,DOWN]-vals[i,j])/abs(vals[i,j])
        print('j_%-8i%-30.16g%-30.16g' % (n,errUP,errDOWN))
        j+=1
    i+=1
    print("\n")
    print("\n")