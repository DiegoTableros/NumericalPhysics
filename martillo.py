#Fisica numerica, tarea 8: Lanzamiento de martillo
#Diego Emilio Dominguez Tableros

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pylab as plt

#Constantes y parametros del modelo
m=7.26                  #Masa de martillo [kg]
rad=0.06                #Radio de martillo [m]
rho=1.2                 #Densidad de aire [kg/m^3]
area=np.pi*rad**2       #Seccion transversal [m^2]
cd_lam=0.5              #Coeficiente rozamiento laminar
cd_inos=0.75            #Coeficiente rozamiento inestable oscilante
coef=[0,cd_lam,cd_inos] #Lista de coeficientes de rozamiento
x0=0                    #Coord x inicial [m]
y0=2                    #Coord y inicial [m]
phi0=np.pi/4            #Angulo de disparo [rad]
dist_rec=86.74          #Distancia record mundial [m]
g=9.81                  #Aceleracion de la gravedad [m/s^2]
#Constantes y parametros de resolucion
h=1e-4          #Tamano de paso para rk4
sol=[]          #Lista/Matriz de soluciones
tiempo=[]       #Lista/Matriz de tiempos
eps=1e-3        #Error absoluto de distancia minimo para calculo de velocidad
ite_max=75      #Maximo de iteraciones para metodo de biseccion
v_aux=32        #Velocidad de prueba auxiliar [m/s]

#Sistema de EDO, dy/dt = f(t,y)
def edo(t,r,cd):
    #x=r[0]; y=r[1]; vx=r[2]; vy=r[3];
    v=np.sqrt(r[2]**2+r[3]**2)          #Rapidez
    f1=r[2]                             #dx/dt=vx
    f2=r[3]                             #dy/dt=vy
    f3=-(rho*area*cd/(2*m))*v*r[2]      #d(vx)/dt=-Fdx
    f4=-g-(rho*area*cd/(2*m))*v*r[3]    #d(vy)/dt=-Fdy-g
    f=np.array((f1,f2,f3,f4))           #Devolver como vector
    return f

#Método de Runge-Kutta de 4to orden para 1 iteracion temporal
def rk4(tn,rn,cd):
    k1=h*edo(tn,rn,cd)
    k2=h*edo(tn+h/2,rn+k1/2,cd)
    k3=h*edo(tn+h/2,rn+k2/2,cd)
    k4=h*edo(tn+h,rn+k3,cd)
    
    r_next=rn+(k1+2*k2+2*k3+k4)/6       #Devolver como vector
    return r_next

#Resolvedor de la EDO con condicion inicial r0 hasta y=0
def solve(r0,cd):
    global sol,tiempo       #Editar las matrices globalmente
    sol=[]; tiempo=[];      #Convertir a listas vacias
    
    r=r0; sol.append([r[0],r[1],r[2],r[3]]);    #1era solucion: r0
    t=0; tiempo.append(t);                      #1ere tiempo: 0 s
    
    #Mientras y>0 m i.e. el martillo no toque el piso (y=0 m)
    while r[1]>0.0:
        t+=h                                #Aumento de tiempo por paso
        r=rk4(t,r,cd)                       #Calculo de nueva solucion
        sol.append([r[0],r[1],r[2],r[3]])   #Adicion a lista de soluciones
        tiempo.append(t)                    #Adicion a lista de tiempos
    
    #Convertir las listas globales a matrices de Numpy
    sol=np.array(sol); tiempo=np.array(tiempo);
    #Retornar la distancia alcanzada en x, alcance horizontal
    return r[0]

#Funcion de determinacion de velocidad para record mundial
def busq_binaria(cd,fig):
    global sol, tiempo      #Editar las matrices globalmente
    
    plt.figure(fig)
    plt.title("Trayectorias calculadas durante la determinación de velocidad inicial para lanzamiento con c_d={}".format(cd))
    print("\nDeterminación de velocidad inicial para lanzamiento con c_d={}\n".format(cd))
    plt.grid(True,which="both"); plt.xlabel("x [m]"); plt.ylabel("y [m]");
    
    v_min=0         #Velocidad minima de busqueda para biseccion [m/s]
    v_max=64        #Velocidad máxima de busqueda para biseccion [m/s]
    ite=1           #Conteo de iteraciones
    
    #Propuesta de velocidad inicial y calculo de alcance horizontal
    v=(v_min+v_max)/2.0
    r0=np.array((x0,y0,v*np.cos(phi0),v*np.sin(phi0)))
    x_final=solve(r0,cd)
    
    #Mientras el alcance horizontal calculado no sea el RECORD MUNDIAl (en <=eps)
    while np.abs(x_final-dist_rec)>eps and ite<ite_max:
        #Impresion en consola del proceso y graficacion
        print("Iteración:",ite,"Alcance horizontal:",x_final,"m")
        plt.plot(sol[:,0],sol[:,1],linewidth=0.5)
        ite+=1
        
        #Metodo de biseccion para acotacion de rango de busqueda de v
        if x_final<=dist_rec:
            v_min=v     #El alcance fue <= que el record, mas velocidad
        else:
            v_max=v     #El alcance fue > que el record, menos velocidad
        
        #Propuesta de velocidad inicial y calculo de alcance horizontal
        v=(v_min+v_max)/2.0
        r0=np.array((x0,y0,v*np.cos(phi0),v*np.sin(phi0)))
        x_final=solve(r0,cd)
    
    #Impresion final para velocidad inicial calculada
    print("Iteración:",ite,"Alcance horizontal:",x_final,"m")
    plt.plot(sol[:,0],sol[:,1],"-k",label="Trayectoria correcta"); plt.legend(); plt.xlim(left=0); plt.ylim(bottom=0);
    print("\nVelocidad inicial calculada con error absoluto de",np.abs(x_final-dist_rec),"al record mundial")
    print(v,"m/s")
    print("------------------------------------------")
    return v

#Determinacion de velocidades iniciales y trayectorias
for i in range(len(coef)):
    #Calculo de velocidad inicial (ya se grafica)
    v_correcta=busq_binaria(coef[i],i+1)
    
    #Grafica de y(t)
    plt.figure(len(coef)+1)
    plt.plot(tiempo,sol[:,1],label="c_d={}".format(coef[i]))
    plt.grid(True,which="both"); plt.legend(); plt.xlabel("t [s]"); plt.ylabel("y [m]"); plt.xlim(left=0); plt.ylim(bottom=0);
    plt.title("Altura del martillo en función del tiempo")
    
    #Grafica de y(x)
    plt.figure(len(coef)+2)
    plt.plot(sol[:,0],sol[:,1],label="c_d={}".format(coef[i]))
    plt.grid(True,which="both"); plt.legend(); plt.xlabel("x [m]"); plt.ylabel("y [m]"); plt.xlim(left=0); plt.ylim(bottom=0);
    plt.title("Trayectoria del martillo")

#Estimacion de la influencia de la distancia en funcion del coef de rozamiento
print("\nEstimación de la influencia del rozamiento en el alcance del tiro para velocidad inicial de",v_aux,"m/s \n")

r0=np.array((x0,y0,v_aux*np.cos(phi0),v_aux*np.sin(phi0)))
coefs=np.arange(0,1.05,0.05)         #Coeficientes de rozamiento
x_frictionless=solve(r0,0)           #El primer alcance sin friccion
alcance=[]                           #El primer cambio en la distancia final es cero

for cd in coefs:
    #Lanzamientos con cond iniciales (x0,y0) y v_aux
    #Obtenemos el alcance horizontal para cada coef de rozamiento
    x_final=solve(r0,cd); alcance.append(x_frictionless-x_final)
    print("Alcance=",x_final,"m  con cd={:.2f}".format(cd))    

#El comportamiento es lineal, grafiquemos el ajuste de datos

def recta(coef,a,b):
    return a*coef+b     #Funcion de ajuste a una recta

print("\nEl cambio en la distancia por rozamiento tiene comportamiento lineal\n")
param,param_cov=curve_fit(recta,coefs,alcance); a=param[0]; b=param[1];
x_aux=np.arange(0,1.01,0.01)
y_aux=a*x_aux+b
print("Parámetros de ajuste para recta: y(x)=ax+b:")
print("a=",a,"b=",b)

#Graficacion del cambio en el alcance en funcion del cd, con Ajuste Lineal
plt.figure(len(coef)+3)
plt.plot(coefs,alcance,'.b',label="(cd,$\Delta$ x)")
plt.plot(x_aux,y_aux,label="Ajuste lineal")
plt.grid(True,which="Both")
plt.legend()
plt.xlabel("Coef de rozamiento (cd)")
plt.ylabel("Cambio en el alcance horizontal ($\Delta$ x) [m]")