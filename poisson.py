#Fisica numerica, tarea 7: Solucion numerica de EDO
#Diego Emilio Dominguez Tableros

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from numpy.linalg import norm
from numpy import inf

x_min=y_min=0               #Limite inferior x,y
x_max=y_max=2*np.pi         #Limite superior x,y
Nx=100                      #Numero de ptos en x
Ny=100                      #Numero de ptos en y
delta_x=(x_max-x_min)/Nx    #Tamaño de paso en x
delta_y=(y_max-y_min)/Ny    #Tamaño de paso en y
delta=delta_x               #Tamaño de paso unificado
tol=3e-2                    #Tolerancia entre normas matriciales
ite_max=150                 #Iteraciones maximas para solucion

#Funcion del lado derecho de la ec. de Poisson
def f(x,y):
    return np.cos(3*x+4*y)-np.cos(5*x-2*y)

#Actualizador de valores de Jacobi
def actualiza_j(k):
    #k es la iteracion del algoritmo
    #k permite guardar la nueva solucion en otra matriz
    for i in range(1,Nx):
        xi=x_min+i*delta_x          #Valor de x actual
        for j in range(1,Ny):
            yi=y_min+j*delta_y      #Valor de y actual
            #Actualizacion de la solucion con sus vecinos
            solJ[i,j,k]=solJ[i+1,j,k-1]+solJ[i-1,j,k-1]+solJ[i,j+1,k-1]+solJ[i,j-1,k-1]
            solJ[i,j,k]-=f(xi,yi)*(delta**2)
            solJ[i,j,k]/=4.0

#Actualizador de valores de Gauss-Seidel
def actualiza_gs(k):
    #Con Gauss-Seidel se trabaja con valores directamente calculados
    solGS[:,:,k]=solGS[:,:,k-1]
    #k es la iteracion del algoritmo
    #k permite guardar la nueva solucion en otra matriz
    for i in range(1,Nx):
        xi=x_min+i*delta_x          #Valor de x actual
        for j in range(1,Ny):
            yi=y_min+j*delta_y      #Valor de y actual
            #Actualizacion de la solucion con sus vecinos
            solGS[i,j,k]=solGS[i+1,j,k]+solGS[i-1,j,k]+solGS[i,j+1,k]+solGS[i,j-1,k]
            solGS[i,j,k]-=f(xi,yi)*(delta**2)
            solGS[i,j,k]/=4.0

#Matriz de mallas para guardar soluciones (Jacobi y Gauss-Seidel)
#Cada malla mide (Nx+1)x(Ny+1)
#Se tiene (ite_max+1) mallas para guardar en cada iteracion
solGS=solJ=np.zeros((Nx+1,Ny+1,ite_max+1),float)
#La condicion inicial para el potencial (malla cero) es aleatoria
solGS[:,:,0]=solJ[:,:,0]=np.random.rand(Nx+1,Ny+1)
##solGS[:,:,0]=solJ[:,:,0]=np.ones((Nx+1,Ny+1))

#Condiciones de frontera
val_ini_x=0.11
val_ini_y=0.24
solJ[:,0,:]=solJ[:,Ny,:]=val_ini_x    #phi(x,0)=phi(x,2pi)
solGS[:,0,:]=solGS[:,Ny,:]=val_ini_x
solJ[0,:,:]=solJ[Nx,:,:]=val_ini_y    #phi(0,y)=phi(2pi,y)
solGS[0,:,:]=solGS[Nx,:,:]=val_ini_y

#Auxiliares para la graficacion de superficie de potencial
X=np.linspace(x_min,x_max,Nx+1)
Y=np.linspace(y_min,y_max,Ny+1)
X,Y=np.meshgrid(Y,X)

#METODO DE JACOBI
############################################################################################
print("Resolucion de la ec. de Poisson por Jacobi")
print("Tolerancia:",tol,"\tNo. max de iteraciones:",ite_max,"\n")
iteJ=1                               #No de iteracion actual (o de malla actual)
norma=norm(solJ[:,:,0],"fro")          #Norma de la cond inicial
actualiza_j(iteJ)                      #Actualizacion de valores
norma_act=norm(solJ[:,:,iteJ],"fro")    #Norma de la primera iteracion
print("Iteración:",iteJ,", diff: ",abs(norma-norma_act))

#Mientras la diferencia entre normas matriciales sea mayor que la tol
while(abs(norma-norma_act)>tol and iteJ<ite_max):
    norma=norm(solJ[:,:,iteJ],"fro")        #Norma de aprox actual
    iteJ+=1                              #Incremento de iteracion
    actualiza_j(iteJ)                      #Actualizacion de valores
    norma_act=norm(solJ[:,:,iteJ],"fro")    #Norma de nueva aprox
    print("Iteración:",iteJ,", diff: ",abs(norma-norma_act))

z_min_J=solJ.min()                     #Valor minimo de potencial
z_max_J=solJ.max()                     #Valor maximo de potencial
figJ=plt.figure(1)
axJ=plt.axes(projection="3d")

#Funcion de graficacion de la cuadro-esima malla solucion
def graficaJ(cuadro,soln,plotJ):
    soln=solJ[:,:,cuadro]        #Malla a graficar
    #Auxiliares de graficacion
    axJ.clear()
    axJ.axes.set_xlim3d(x_min,x_max) 
    axJ.axes.set_ylim3d(y_min,y_max) 
    axJ.axes.set_zlim3d(z_min_J,z_max_J)
    plotJ=axJ.plot_surface(X,Y,soln,cmap="plasma",antialiased=False,linewidth=0)
    axJ.set_xlabel('x')
    axJ.set_ylabel('y')
    axJ.set_zlabel('$\phi$(x,y)')
    axJ.axes.set_title(f"Potencial $\phi(x,y)$ por Jacobi.\nCuadro: {cuadro}")
    return plotJ

#Primer cuadro de animacion y generacion de animacion
plotJ=axJ.plot_surface(X,Y,solJ[:,:,0],cmap="plasma",antialiased=False,linewidth=0)
aniJ=animation.FuncAnimation(figJ,graficaJ,frames=iteJ,fargs=(solJ,plotJ),interval=0,blit=False)

#METODO DE GAUSS-SEIDEL
############################################################################################
print("\n\nResolucion de la ec. de Poisson por Gauss-Seidel")
print("Tolerancia:",tol,"\tNo. max de iteraciones:",ite_max,"\n")
iteGS=1                               #No de iteracion actual (o de malla actual)
norma=norm(solGS[:,:,0],"fro")          #Norma de la cond inicial
actualiza_gs(iteGS)                      #Actualizacion de valores
norma_act=norm(solGS[:,:,iteGS],"fro")    #Norma de la primera iteracion
print("Iteración:",iteGS,", diff: ",abs(norma-norma_act))

#Mientras la diferencia entre normas matriciales sea mayor que la tol
while(abs(norma-norma_act)>tol and iteGS<ite_max):
    norma=norm(solGS[:,:,iteGS],"fro")        #Norma de aprox actual
    iteGS+=1                              #Incremento de iteracion
    actualiza_gs(iteGS)                      #Actualizacion de valores
    norma_act=norm(solGS[:,:,iteGS],"fro")    #Norma de nueva aprox
    print("Iteración:",iteGS,", diff: ",abs(norma-norma_act))

z_min_GS=solGS.min()                     #Valor minimo de potencial
z_max_GS=solGS.max()                     #Valor maximo de potencial
figGS=plt.figure(2)
axGS=plt.axes(projection="3d")

#Funcion de graficacion de la cuadro-esima malla solucion
def graficaGS(cuadro,soln,plotGS):
    soln=solGS[:,:,cuadro]        #Malla a graficar
    #Auxiliares de graficacion
    axGS.clear()
    axGS.axes.set_xlim3d(x_min,x_max) 
    axGS.axes.set_ylim3d(y_min,y_max) 
    axGS.axes.set_zlim3d(z_min_GS,z_max_GS)
    plotGS=axGS.plot_surface(X,Y,soln,cmap="viridis",antialiased=False,linewidth=0)
    axGS.set_xlabel('x')
    axGS.set_ylabel('y')
    axGS.set_zlabel('$\phi$(x,y)')
    axGS.axes.set_title(f"Potencial $\phi(x,y)$ por Gauss-Seidel.\nCuadro: {cuadro}")
    return plotGS

#Primer cuadro de animacion y generacion de animacion
plotGS=axGS.plot_surface(X,Y,solGS[:,:,0],cmap="viridis",antialiased=False,linewidth=0)
aniGS=animation.FuncAnimation(figGS,graficaGS,frames=iteGS,fargs=(solGS,plotGS),interval=0,blit=False)