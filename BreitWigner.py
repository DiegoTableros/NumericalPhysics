#Tarea 6, Fisica Numerica. Diego Emilio Dominguez Tableros
#Problema 3. Ajuste de Breit-Wigner

import numpy as np
import matplotlib.pyplot as plt

#Lectura de datos de archivo BretWigner.txt
datos=np.loadtxt("BretWigner.txt",skiprows=1,dtype=float)
E=datos[:,1]        #Datos de energia E [MeV]
f_E=datos[:,2]      #Datos de f(E) [MeV]
sigma=datos[:,3]    #Incertidumbres de f(E)
n_datos=len(E)      #Numero de datos
h=1e-4              #Paso para aproximaciones
tolerancia=1e-6     #Tolerancia absoluta entre aproximaciones

#Funciones para Newton-Raphson, pasar a linea 105
######################################################################
def Derivada(f,vec):
    #Aproximación de la derivada de una función por forward difference
    #Puntos de evaluacion de la derivada x, y, z
    x=vec[0]; y=vec[1]; z=vec[2]
    #Aproximacion de derivadas parciales en cada variable
    dfx=(f(x+h,y,z)-f(x,y,z))/h
    dfy=(f(x,y+h,z)-f(x,y,z))/h
    dfz=(f(x,y,z+h)-f(x,y,z))/h
    #Construccion de vector de salida: f'(vec)
    df=np.array([dfx,dfy,dfz]); df.shape=(3,1);
    return df

def MatrizJacobiana(f1,f2,f3,vec):
    #Matriz jacobiana de las funciones f1, f2, f3 evaluada en vec
    df1=Derivada(f1,vec)
    df2=Derivada(f2,vec)
    df3=Derivada(f3,vec)
    #Construccion de matriz Jacobiana de salida: J(vec)
    J=np.array([df1,df2,df3]); J.shape=(3,3);
    return J

def evaluar_f(f1,f2,f3,vec):
    #Evaluacion de la funcion f=[f1(x,y,z),f2(x,y,z),f3(x,y,z)] en vec
    f1_eval=f1(vec[0],vec[1],vec[2])
    f2_eval=f2(vec[0],vec[1],vec[2])
    f3_eval=f3(vec[0],vec[1],vec[2])
    #Construccion de vector de salida: f(vec)
    f_eval=np.array([f1_eval,f2_eval,f3_eval]); f_eval.shape=(3,1);
    return f_eval    

def NewtonRaphsonMulti(f1,f2,f3,semilla,tol):
    #Implementacion de metodo NewtonRaphson multivariable
    #f1,f2,f3 son las funciones cuyas raices se quieren encontrar
    #es decir, se desea encontrar f=[f1,f2,f3]=[0,0,0]=0
    #semilla es un vector con valores cercanos a dichas raices
    #La norma vectorial usada es la de Frobenius (aka. Euclidea 2)
    print("Aproximacion de raices por Newton-Raphson con tolerancia de ",tol)
    
    #Primera iteracion para aproximar raiz
    ite=0           #Conteo de iteraciones
    r=semilla
    print("x_",ite,"=\n",r,"\n"); ite+=1;
    J=MatrizJacobiana(f1,f2,f3,r)
    J_inversa=np.linalg.inv(J)
    f_eval=evaluar_f(f1,f2,f3,r)
    raiz=r-np.matmul(J_inversa,f_eval)
    
    #Mientras que las aprox difieran en mas de la tolerancia
    #se debe seguir calculando otra aproximacion
    while(np.linalg.norm(raiz-r)>tol):
        r=raiz                                  #Actualiza valor de aprox
        print("x_",ite,"=\n",r,"\n"); ite+=1;   #Impresion de iteracion
        J=MatrizJacobiana(f1,f2,f3,r)           #Matriz jacobiana evaluada
        J_inversa=np.linalg.inv(J)              #Inversa de matriz jacobiana
        f_eval=evaluar_f(f1,f2,f3,r)            #Funcion evaluada
        raiz=r-np.matmul(J_inversa,f_eval)      #x_sig = x_act - (J_inv * f_eval)
    
    print("x_",ite,"=\n",r,"\n");               #Impresion final de raiz
    return raiz

def g(a1,a2,a3,x):
    #Evaluacion de g(x) donde x puede pertenecer a los datos E
    return a1/((x-a2)**2+a3)
  
def f_1(a1,a2,a3):
    #Primera funcion para ajuste por minimizacion de chi cuadrada
    valor=0
    for i in range(n_datos):
        x=E[i]; y=f_E[i]; s=sigma[i];   #Valores para termino
        valor+=(y-g(a1,a2,a3,x))/(((x-a2)**2+a3)*(s**2))
    return valor

def f_2(a1,a2,a3):
    #Segunda funcion para ajuste por minimizacion de chi cuadrada
    valor=0
    for i in range(n_datos):
        x=E[i]; y=f_E[i]; s=sigma[i];   #Valores para termino
        valor+=(y-g(a1,a2,a3,x))*(x-a2)/((((x-a2)**2+a3)**2)*(s**2))
    return valor

def f_3(a1,a2,a3):
    #Segunda funcion para ajuste por minimizacion de chi cuadrada
    valor=0
    for i in range(n_datos):
        x=E[i]; y=f_E[i]; s=sigma[i];   #Valores para termino
        valor+=(y-g(a1,a2,a3,x))/((((x-a2)**2+a3)**2)*(s**2))
    return valor
#############################################################

#Cálculo de raíces para ajuste por NewtonRaphson multivariable

#La semilla inicial debe estar cerca de los valores esperados
seed_inicial=np.array([60000,50,700])
seed_inicial.shape=(3,1)

raices=NewtonRaphsonMulti(f_1,f_2,f_3,seed_inicial,tolerancia)
a1=raices[0]            #Coeficiente de ajuste a1
a2=raices[1]            #Coeficiente de ajuste a2
a3=raices[2]            #Coeficiente de ajuste a3
print("Parametros de ajuste para Breit-Weigner")
print("a1=",a1,"\na2=",a2,"\na3=",a3)

#Graficacion
x_aux=np.linspace(min(E),max(E),250)
y_aux=g(a1,a2,a3,x_aux)
plt.plot(x_aux,y_aux,label="Ajuste Min. Cuadrados")
#Graficacion de datos originales con barras de error para f(E)
error_f_E=[sigma,sigma]
plt.errorbar(E,f_E,yerr=error_f_E,fmt='.',label="Datos")

plt.xlabel("Energía E [MeV]")
plt.ylabel("f(E) [MeV]")
plt.legend()
plt.grid(which="both")
plt.show()