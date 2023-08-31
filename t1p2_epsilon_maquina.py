#Tarea 1, problema 2. Presición de la máquina
#Diego Emilio Dominguez Tableros

#El eps_m se define como la cantidad más grande
#para la cual 1.0+eps_m=1.0

#Partimos de eps_m=1
eps_m=1.0
#Encabezado de impresión
print("%-24s%-24s" % ("eps_m","1.0 + eps_m"))
print("------------------------------------------------")

#Mientras 1.0+eps_m sea diferente de 1.0
while (1.0+eps_m)!=1.0:
    #Dividimos el valor de eps_m a la mitad
    eps_m/=2
    #Imprimimos el eps_m actual
    print("%-24.16g%-24.16f" % (eps_m,1.0+eps_m))

#Al finalizar, eps_m contiene el epsilon de la máquina dentro
#de un factor de 2 como consecuencia de ir diviendo entre 2
#en cada iteración

print("\nEpsilón de máquina (eps_m): ",eps_m)