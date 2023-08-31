#Estudio de la evolucion de la temperatura en una barra

from numpy import *
import matplotlib.pylab as p
from mpl_toolkits.mplot3d import Axes3D

#Constantes
L=1         #Longitud de la barra [m]
t_f=4000    #Tiempo final o total [s]
Nx=100      #Numero de pasos para x
Nt=10000    #Numero de pasos para t
dx=L/Nx     #Tamano de paso para x
dt=t_f/Nt   #Tamano de paso para t
T_0=100     #Temperatura inicial de barra [∞C]
T_f=0       #Temperatura en los extremos de barra [∞C]
rho=2698    #Densidad de barra (Al) [kg/m^3]
K=237       #Conductividad termica (Al) [W/(∞C*m)]
C=897       #Calor especifico (Al) [J/(C*kg)]
kappa=K/(C*rho)         #Coeficiente de difusion
eta=kappa*dt/(dx**2)    #Factor eta

print("Factor eta: ", eta)  #Condicion estabilidad VonNeunmman-Courant
print("Tama√±o de paso dx: ", dx);
print("Tama√±o de paso dt: ", dt);
T=zeros((Nt+1,Nx+1),float)  #Matriz de soluciones (malla)

#Inicializacion de condiciones iniciales/de frontera
#Condicion de frontera T(x=0,t)=T(x=L,t)=T_f, 0<=t<=t_f
for j in range(0,Nt+1):
    T[j,0]=T_f
    T[j,Nx]=T_f    
#Condicion inicial T(x,t=0)=T_0, 0<x<L
for i in range(1,Nx):
    T[0,i]=T_0

#Resolucion de la EDP sobre la malla con dif. finitas
for j in range(0,Nt):
    for i in range (1,Nx):
        T[j+1,i]=T[j,i]+eta*(T[j,i+1]+T[j,i-1]-2*T[j,i])

#Graficacion de solucion numerica: FIGURA 1
x_aux=list(range(0,Nx+1))
y_aux=list(range(0,Nt+1))
Xg,Yg=p.meshgrid(x_aux,y_aux)
Zg=T
fig=p.figure()
ax=Axes3D(fig)
ax.plot_wireframe(Yg,Xg,Zg)
ax.set_ylabel("Posicion x")
ax.set_xlabel("Tiempo t")
ax.set_zlabel("Temperatura T(x,t)")
p.title("Evolucion de la temperatura de una barra")
p.show()

#Solucion analitica
n_term_max=200      #Numero max de terminos de serie
x_aux=p.arange(0,L,L/100)
y_aux=p.arange(0,t_f,t_f/1000)
Xg,Yg=p.meshgrid(x_aux,y_aux)
Zg=zeros((len(y_aux),len(x_aux)),float)
#Funcion calculadora de terminos de serie
def term_sol_ana(n):
    term=200
    term*=(1-cos(n*pi))/(n*pi)
    term*=sin(n*pi*Xg)
    term*=exp(-(n**2)*(pi**2)*kappa*Yg)
    return term
#Calculo de solucion analitica
for n in range(1,n_term_max+1):
    term=term_sol_ana(n)
    Zg=Zg+term
#Graficacion de solucion analitica: FIGURA 2
fig2=p.figure()
ax=Axes3D(fig2)
ax.plot_wireframe(Yg,Xg,Zg)
ax.set_ylabel("Posicion x [m]")
ax.set_xlabel("Tiempo t [s]")
ax.set_zlabel("Temperatura T(x,t) [∞C]")
p.title("Solucion analitica de T(x,t)")
p.show()