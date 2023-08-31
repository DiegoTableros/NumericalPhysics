#Estudio de las vibraciones en una cuerda

from numpy import *
import matplotlib.pylab as p
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.pyplot as plt

#Constantes
L=1         #Longitud de la cuerda [m]
t_f=5       #Tiempo final o total [s]
Nx=50      #Numero de pasos para x
Nt=3000    #Numero de pasos para t
dx=L/Nx     #Tamano de paso para x
dt=t_f/Nt   #Tamano de paso para t
y0=0.001      #Amplitud inicial en x=x0 [m]
x0=0.3      #Pos. de amplitud inicial triangular [m]
            #WARNING: x0 debe ser multiplo entero de dx
Amp_0=0     #Amplitud en los extremos [m]
T=10        #Tension en la cuerda [N]
rho=0.1     #Densidad lin. de masa [kg/m]
c=sqrt(T/rho)       #Vel. de propagacion [m/s]
c_prima=dx/dt       #Vel. de la malla          
p_estable=c/c_prima #Parametro de estabilidad Courant

print("Parametro de estabilidad: ", p_estable)  #Condicion estabilidad Courant
print("Tamano de paso dx: ", dx);
print("Tamano de paso dt: ", dt);
y=zeros((Nt+1,Nx+1),float)  #Matriz de soluciones (malla)

#Inicializacion de condiciones iniciales/de frontera
#Condicion de frontera y(x=0,t)=y(x=L,t)=Amp_0, 0<=t<=t_f
for j in range(0,Nt+1):
    y[j,0]=Amp_0
    y[j,Nx]=Amp_0
#Condicion inicial y(x,t=0)=y0, 0<x<L
i0=int(x0/dx)          #Indice corresp. a x=x0
for i in range(1,Nx):
    if i<=i0:
        y[0,i]=(y0/x0)*dx*i             #x<=x0
    else:
        y[0,i]=(y0/(x0-L))*(dx*i-L)     #x>x0
#Condicion inicial dy/dt=0 para t=dt
for i in range(1,Nx):
    y[1,i]=y[0,i]       #La solución es la anterior

#Resolucion de la EDP sobre la malla con dif. finitas
fac=(c**2)/(c_prima**2)
for j in range(1,Nt):
    for i in range (1,Nx):
        #Resolucion por diferencias finitas
        y[j+1,i]=2*y[j,i]-y[j-1,i]+fac*(y[j,i+1]+y[j,i-1]-2*y[j,i])

#Graficacion
fig, ax = plt.subplots()
x=arange(0,L+dx,dx)
equis=list(range(0,Nx+1))   #Rescatador de valores
line,=ax.plot(x,y[0,equis])

def animate(i):
    line.set_ydata(y[i,equis])  #Actualización de datos y
    return line,

ani = animation.FuncAnimation(fig,animate,interval=1,blit=False)
plt.xlim([0,L])
plt.ylim([-abs(2*y0),abs(2*y0)])
ax.set_ylabel("Amplitud y(x,t) [m]")
ax.set_xlabel("Eje x [m]")
plt.title("Vibraciones en cuerda con extremos fijos")
plt.grid(visible=True,which="both")
plt.show()