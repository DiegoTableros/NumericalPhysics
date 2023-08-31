#Tarea 1, problema 1. Overflow y Underflow
#Diego Emilio Dominguez Tableros

#Para determinar el límite del Underflow, se debe encontrar el
#número más pequeño que sea considerado cero por Python

#Partiendo de +1
aux_under=1.0
#Mientras este número sea diferente de cero
while aux_under!=0.0:
    #Se divide el número entre 2 y se guarda el actual
    under=aux_under
    aux_under/=2.0
#Al salir del ciclo, under contiene el límite del bajo flujo dentro de un factor de 2
print("Límite de underflow positivo:")
print("%-24.16g" % (under))

#Partiendo de -1
aux_under=-1.0
#Mientras este número sea diferente de cero
while aux_under!=0.0:
    #Se divide el número entre 2 y se guarda el actual
    under=aux_under
    aux_under/=2.0
#Al salir del ciclo, under contiene el límite del bajo flujo dentro de un factor de 2
print("Límite de underflow negativo:")
print("%-24.16g" % (under))

print("\n");

#Para determinar el límite del Overflow, se debe encontrar el
#número más grande que Python no considere infinito

#Partiendo de +1
aux_over=1.0
#Mientras este numero
while aux_over!=float("inf"):
    #Se multiplica el número por 2 y se guarda el actual
    over=aux_over
    aux_over*=2
    #print(over,"\n");
#Al salir del ciclo, over contiene el limite del sobre flujo dentro de un factor de 2
print("Límite de overflow positivo:")
print("%-24.16g" % (over))

#Partiendo de -1
aux_over=-1.0
#Mientras este numero
while aux_over!=-float("inf"):
    #Se multiplica el número por 2 y se guarda el actual
    over=aux_over
    aux_over*=2
    #print(over,"\n");
#Al salir del ciclo, over contiene el limite del sobre flujo dentro de un factor de 2
print("Límite de overflow negativo:")
print("%-24.16g" % (over))