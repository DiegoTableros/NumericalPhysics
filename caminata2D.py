#Tarea 5. Fisica Numerica. Caminatas aleatorias en 2D
#Diego Emilio Dominguez Tableros

import numpy as np
import matplotlib.pylab as plt

def caminata_2D(n):
    k=int(np.sqrt(n))               #Numero de experimentos
    dist_cua=np.zeros(k)            #Arreglo de R^2 para c/experimento
    for expe in range(k):
        #Nuevo experimento "expe"
        tray=np.zeros((n+1,2),float)            #Arreglo de trayectoria
        pasos=np.random.uniform(-1,1,(n,2))     #Arreglo de n pasos
        for i in range(n):
            L=np.sqrt(pasos[i,0]**2+pasos[i,1]**2)  #Longitud de paso
            pasos[i,0]/=L                           #Normalizacion X
            pasos[i,1]/=L                           #Normalizacion Y
            tray[i+1,0]=tray[i,0]+pasos[i,0]        #Actualiza tray X
            tray[i+1,1]=tray[i,1]+pasos[i,1]        #Actualiza tray Y
            
        dist_cua[expe]=tray[n,0]**2+tray[n,1]**2    #R^2 para este expe
    
    prom_dist_cua=sum(dist_cua)/k           #Promedio R^2 para K expe
    R_rms=np.sqrt(prom_dist_cua)            #RootMeanSquare para este exp
    
    #Hipotesis teorica ( <x_ii * x_jj>/R^2 ) approx 0
    ii=1
    jj=np.random.randint(2,n)
    print("Hipotesis teorica para N=",n,"pasos:")
    print(((pasos[ii,0]*pasos[jj,0])/n)/prom_dist_cua**2)
    
    plt.figure(1)
    plt.plot(tray[:,0],tray[:,1],"-",linewidth=0.5)         #Grafica trayectoria
    
    return R_rms                            #Retorna R_rms de este N

N=1000                                        #N limite para experimentos (grande)
paso=int(N/10)                                  #Numero de pruebas de caminata
plt.figure(2)
plt.plot([1,np.sqrt(N)],[1,np.sqrt(N)],"-")     #Funcion identidad
for i in range(paso,N+1,paso):
    #Obtenemos R_rms de esta caminata
    R_rms_i=caminata_2D(i)
    #Graficamos (sqrt(N),R_rms(N)) para esta caminata
    plt.figure(2)
    plt.plot(np.sqrt(i),R_rms_i,'.')

#Graficacion
plt.figure(1)
plt.xlabel("x")
plt.ylabel("y")
plt.plot(0,0,'.',color="black",markersize=10)
plt.figure(2)
plt.xlabel("sqrt(N)")
plt.ylabel("R_rms(N)")
plt.show()