#Tarea 6, Fisica Numerica. Diego Emilio Dominguez Tableros
#Problema 5. Espectro cuerpo negro y COBE

import numpy as np
import matplotlib.pyplot as plt

#Lectura de datos de archivo COBE.txt
datos=np.loadtxt("COBE.txt",skiprows=1,dtype=float)
fre=datos[:,0]      #Datos de frecuencia [1/cm]
pot=datos[:,1]      #Datos de potencia [MJy/sr]
sigma=datos[:,2]    #Incertidumbres
sigma/=1000         #Escalamiento de unidades
n_datos=len(fre)    #Numero de datos
h=1e-5              #Paso para aproximaciones
tolerancia=1e-7     #Tolerancia absoluta entre aproximaciones

#Funciones para Newton-Raphson, pasar a linea 105
######################################################################
def Derivada(f,vec):
    #Aproximación de la derivada de una función por forward difference
    #Puntos de evaluacion de la derivada x, y
    x=vec[0]; y=vec[1];
    #Aproximacion de derivadas parciales en cada variable
    dfx=(f(x+h,y)-f(x,y))/h
    dfy=(f(x,y+h)-f(x,y))/h
    #Construccion de vector de salida: f'(vec)
    df=np.array([dfx,dfy]); df.shape=(2,1);
    return df

def MatrizJacobiana(f1,f2,vec):
    #Matriz jacobiana de las funciones f1, f2, f3 evaluada en vec
    df1=Derivada(f1,vec)
    df2=Derivada(f2,vec)
    #Construccion de matriz Jacobiana de salida: J(vec)
    J=np.array([df1,df2]); J.shape=(2,2);
    return J

def evaluar_f(f1,f2,vec):
    #Evaluacion de la funcion f=[f1,f2] en vec
    f1_eval=f1(vec[0],vec[1])
    f2_eval=f2(vec[0],vec[1])
    #Construccion de vector de salida: f(vec)
    f_eval=np.array([f1_eval,f2_eval]); f_eval.shape=(2,1);
    return f_eval    

def NewtonRaphsonMulti(f1,f2,semilla,tol):
    #Implementacion de metodo NewtonRaphson multivariable
    #f1,f2,f3 son las funciones cuyas raices se quieren encontrar
    #es decir, se desea encontrar f=[f1,f2]=[0,0]=0
    #semilla es un vector con valores cercanos a dichas raices
    #La norma vectorial usada es la de Frobenius (aka. Euclidea 2)
    print("Aproximacion de raices por Newton-Raphson con tolerancia de ",tol)
    
    #Primera iteracion para aproximar raiz
    ite=0           #Conteo de iteraciones
    r=semilla
    print("x_",ite,"=\n",r,"\n"); ite+=1;
    J=MatrizJacobiana(f1,f2,r)
    J_inversa=np.linalg.inv(J)
    f_eval=evaluar_f(f1,f2,r)
    raiz=r-np.matmul(J_inversa,f_eval)
    
    #Mientras que las aprox difieran en mas de la tolerancia
    #se debe seguir calculando otra aproximacion
    while(np.linalg.norm(raiz-r)>tol):
        r=raiz                                  #Actualiza valor de aprox
        print("x_",ite,"=\n",r,"\n"); ite+=1;   #Impresion de iteracion
        J=MatrizJacobiana(f1,f2,r)           #Matriz jacobiana evaluada
        J_inversa=np.linalg.inv(J)              #Inversa de matriz jacobiana
        f_eval=evaluar_f(f1,f2,r)            #Funcion evaluada
        raiz=r-np.matmul(J_inversa,f_eval)      #x_sig = x_act - (J_inv * f_eval)
    
    print("x_",ite,"=\n",r,"\n");               #Impresion final de raiz
    return raiz

def g(a1,a2,x):
    #Evaluacion de g(x) donde x puede pertenecer a los datos E
    return (a1*(x**3))/(np.exp(a2*x)-1)
  
def f_1(a1,a2):
    #Primera funcion para ajuste por minimizacion de chi cuadrada
    valor=0
    for i in range(n_datos):
        x=fre[i]; y=pot[i]; s=sigma[i];   #Valores para termino
        valor+=((y-g(a1,a2,x))/s**2)*(x**3/(np.exp(a2*x)-1))
    return valor

def f_2(a1,a2):
    #Segunda funcion para ajuste por minimizacion de chi cuadrada
    valor=0
    for i in range(n_datos):
        x=fre[i]; y=pot[i]; s=sigma[i];   #Valores para termino
        valor+=((y-g(a1,a2,x))/s**2)*((-a1*a2*(x**3)*np.exp(a2*2))/((np.exp(a2*x)-1)**2))
    return valor

#############################################################

#Graficacion de datos originales con barras de error
error_volt=[sigma,sigma]
plt.errorbar(fre,pot,yerr=error_volt,fmt='.',label="Datos")

#Cálculo de raíces para ajuste por NewtonRaphson multivariable

#La semilla inicial debe estar cerca de los valores esperados
seed_inicial=np.array([34,0.5])
seed_inicial.shape=(2,1)

raices=NewtonRaphsonMulti(f_1,f_2,seed_inicial,tolerancia)
a1=raices[0]            #Coeficiente de ajuste a1
a2=raices[1]            #Coeficiente de ajuste a2
print("Parametros de ajuste para decaimiento de voltaje")
print("a1=",a1,"\na2=",a2,"\n")

#Graficacion: Potencia en funcion de frecuencia
x_aux=np.linspace(min(fre),max(fre),500)
y_aux=g(a1,a2,x_aux)
plt.plot(x_aux,y_aux,label="Ajuste Min. Cuadrados")

plt.xlabel("Frecuencia [1/cm]")
plt.ylabel("Potencia radiada [MJy/sr]")
plt.legend()
plt.grid()
plt.show()

#Calculo de coeficiente Chi cuadrada
chi_cua=0
for i in range(n_datos):
    chi_cua+=((pot[i]-g(a1,a2,fre[i]))/sigma[i])**2
print("Coeficiente Chi cuadrada: ",chi_cua)
print("Numero de datos: ",n_datos)