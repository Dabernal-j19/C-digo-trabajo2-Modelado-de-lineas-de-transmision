#Importacion de librerias
import numpy as np #para todo
import matplotlib.pyplot as plt   #para las graficas
#Importacion de librerias




#-----------------------
#DEFINICION DE FUNCIONES
#-----------------------
#Creacion de la "matriz de distancias"
def calcular_distancias(coordenadas_conductores,cantidad_conductores):  #recordar que se se asume que los conductores son de radios aproximadamente iguales, y que la distancia entre conductores es mucho mayor que el radio de los mismos. Es por ello que se calcula la distancia con el centro del conductor
    coordenadas = np.array(coordenadas_conductores)
    distancias = np.zeros((cantidad_conductores, cantidad_conductores))
    for i in range(cantidad_conductores):
        for j in range(cantidad_conductores):
            if i != j:
                distancias[i, j] = np.linalg.norm(coordenadas[i] - coordenadas[j])
    return distancias
#Creacion de la "matriz de distancias"


#Resistencia general (con o sin efecto piel)
def resistencia_general(radio, resistividad_conductor, delta):
    if delta > radio:  #Si el efecto piel es despreciable (delta > radio), usar la fórmula clásica
        mensaje1="Por la frecuencia de la linea, delta es mayor que el radio de los conductores. Debido a lo anterior NO se tiene en cuenta el efecto piel"
        area_transversal = np.pi * radio**2
        resistencia_metro = resistividad_conductor / area_transversal
    else:  #Si el efecto piel es dominante (delta <= radio), usar la fórmula con efecto piel
        mensaje1="Por la frecuencia de la linea, delta es menor o igual que el radio de los conductores. Debido a lo anterior SI se tiene en cuenta el efecto piel"
        resistencia_metro = resistividad_conductor / (2 * np.pi * radio * delta)
    return resistencia_metro, mensaje1
#Resistencia general (con o sin efecto piel)


#Inductancia propia apriximada (interna teniendo en cuenta el efecto piel y teniendo en cuenta el retorno por tierra)
def inductancia_propia_aprox_efecto_piel_tierra_y_sin_efecto_piel(r,mu_medio,altura_conductor,p,mu_conductor,resistividad_conductor,frecuencia):
    if p==0:    #Este condicional es para saber si el terreno tiene conductividad infinita, calcular la inductancia normalmente
        inductancia_final=(mu_conductor/8*np.pi)+(mu_medio*np.log(((2*altura_conductor)/r)))/(2*np.pi)  #formula para la inductancia normal mas la inductancia interna sin tener en cuenta el efecto piel
    else:
        inductancia_interna_efecto_piel=(1/(4*np.pi*r))*(np.sqrt(mu_conductor/(np.pi*(1/resistividad_conductor))))*(np.sqrt(1/frecuencia))
        inductancia_externa_retorno_tierra=(mu_medio/(2*np.pi))*np.log((altura_conductor+p)/r)
        inductancia_final=inductancia_externa_retorno_tierra+inductancia_interna_efecto_piel
    return inductancia_final
#Inductancia propia apriximada (interna teniendo en cuenta el efecto piel y teniendo en cuenta el retorno por tierra)


#Inductancia mutua aproximada
def inductancia_mutua_aproximada(mu_medio,distancia_entre_conductores,distancia_imagen):
    inductancia_mutua=(mu_medio*np.log(distancia_imagen/distancia_entre_conductores))/(2*np.pi)   #calcula la inductancia mutua entre dos conductores
    return inductancia_mutua
#Inductancia mutua aproximada


#Funcion para graficar
def graficar2ConjuntosDeDatos(t,conjuntoDatos1,nombre1,conjuntoDatos2,nombre2,tituloX,tituloY,tituloGrafica,forma1,forma2):
  plt.figure(figsize=(8, 5))
  plt.plot(t,conjuntoDatos1,forma1,color="green",label=nombre1)
  plt.plot(t,conjuntoDatos2,forma2,color="black",label=nombre2)
  plt.xlabel(tituloX)
  plt.ylabel(tituloY)
  plt.title(f"{tituloGrafica}")
  plt.grid(True)
  plt.legend()
  plt.tight_layout()
  plt.show()
#Funcion para graficar


#Funcion para generar la matriz de Forstcue segun la cantidad de circuitos trifasicos
def generar_matriz_fortescue(cantidad_circuitos_trifasicos):
    a = (-1/2) + (np.sqrt(3)/2)*1j  #Definición del número complejo a=1|120°
    a2 = a**2
    T_base = np.array([[1, 1, 1],    
                       [1, a2, a],
                       [1, a, a2]])   #Matriz base 3x3 (para un circuito trifásico)
    T = np.zeros((3*cantidad_circuitos_trifasicos, 3*cantidad_circuitos_trifasicos), dtype=complex)   #Construcción de la matriz bloque-diagonal para n circuitos trifásicos
    for i in range(cantidad_circuitos_trifasicos):
        inicio = 3 * i
        fin = inicio + 3
        T[inicio:fin, inicio:fin] = T_base
    T_inv = np.linalg.inv(T)   #Calcular inversa
    return T, T_inv
#Funcion para generar la matriz de Forstcue segun la cantidad de circuitos trifasicos
#-----------------------
#DEFINICION DE FUNCIONES
#-----------------------




#------------------------
#DIGITACION DE PARAMETROS
#------------------------
#Coordenadas de los conductores los valores x e y deben estar en metros [m]. Los ultimos valores deben ser los cables de guarda. Tener en cuenta como es el plano cartesiano de referencia que se propuso
coordenadas_conductores=np.array([(0,20),      
                                  (0.6,22.5),
                                  (0,25),
                                  (0,28)])

"""coordenadas_conductores=np.array([(2.5,20),
                                  (3.1,22.5),
                                  (2.5,25),
                                  (-2.5,20),
                                  (-3.1,22.5),
                                  (-2.5,25),
                                  (2.5,28),
                                  (-2.5,28)]) """
#Coordenadas de los conductores los valores x e y deben estar en metros [m]. Los ultimos valores deben ser los cables de guarda. Tener en cuenta como es el plano cartesiano de referencia que se propuso


#Otros parametros
cantidad_cables_guarda=1  #candidad de conductores que toman el papel de cables de guarda
cantidad_circuitos_trifasicos=1   #cantidad de circuitos trifasicos de la linea de transmision. Este programa SOLO PUEDE TRABAJAR CON  n CIRCUITOS TRIFASICOS. Esto se usa pra generar la matriz de Fortscue
voltaje_linea=115000  #valor en voltios [V]
radio_conductores_fases=0.00598  #valor en metros [m]
radio_cable_guarda=0.0041   #valor en metros [m]
longitud_linea=1000   #valor en metros [m]
epsilon_r_medio=1   #Permitividad relativa del medio circundante de la linea de transmision (epsilon_r) [adimensional]
mu_r_medio=1  #Permeabilidad relativa del medio circundante de la linea de transmision (mu_r) [adimensional]
sigma_medio=3e-8 #Conductividad del medio circundante de la linea de transmision [S/m]
sigma_tierra="inf" #Conductividad de la tierra [S/m]. Este valor es para tener en cuenta el efecto de penetracion que se genera por usarla de retorno. En caso que no se quiera tener en cuenta este efecto (tener conductividad infinita), colocar textualmente: "inf"
mu_r_conductor=1.000022 #Permeabilidad relativa del conductor de la linea de transmision (mu_r) [adimensional]
resistividad_conductor=2.82e-8   #Resistividad de los conductores [Ohm*m]
frecuencia=60  #Frecuencia de la fuente de voltaje [Hz]
#Otros parametros
#------------------------
#DIGITACION DE PARAMETROS
#------------------------




#---------------------------------------------------------------
#CALCULOS PREVIOS Y GRAFICA DEL SISTEMA DE CONDUCTORES INGRESADO
#---------------------------------------------------------------
#Calculos preelimnares y declaracion de constantes necesarias
epsilon_0=8.854e-12 #epsilon_0=8.854e-12 [F/m]
mu_0=(np.pi)*4e-7 #mu_0=1.2566e-6 [H/m]
epsilon_medio=epsilon_r_medio*epsilon_0   #epsilon medio absoluto
mu_medio=mu_r_medio*mu_0   #mu medio absoluto
mu_conductor=mu_r_conductor*mu_0  #mu conductor absoluto
w=2*np.pi*frecuencia   #frecuencia angular
delta=np.sqrt((2*resistividad_conductor)/(w*mu_conductor))  #Delta, constante para tener en cuenta el "skin depth" o el efecto piel
np.set_printoptions(linewidth=400)  #Esto ajusta en la terminal para poder mostrar toda la fila de la matriz en una sola columna

if sigma_tierra!="inf":   #Este condicional es para saber si tener en cuenta o no el efecto de penetracion en la tierra
    p=1/np.sqrt(1j*w*mu_medio*sigma_tierra)  #p, distancia compleja para calcular la inductancia generada por utilizar la tierra como retorno
else:
    p=0

cantidad_conductores=coordenadas_conductores.shape[0] #cantidad de conductores
#Calculos preelimnares y declaracion de constantes necesarias


#Extraer coordenadas x e y de cada conductor, y generar un vector para la tierra apartir de las coordenadas en x
coordenadasX=coordenadas_conductores[:,0]  #Array para almacenar todas las coordenadas x
coordenadasY=coordenadas_conductores[:,1]  #Array para almacenar todas las coordenadas y
gnd=np.zeros(cantidad_conductores)  #genera un array del mismo tamaño que coordenadasX, con todos sus valores iguales a 0
#Extraer coordenadas x e y de cada conductor, y generar un vector para la tierra apartir de las coordenadas en x


#Calculo de distancias entre conductores
matriz_distancias=calcular_distancias(coordenadas_conductores,cantidad_conductores)
#Calculo de distancias entre conductores


#Introducir la altura de cada conductor en la diagonal de la matriz de distancias. Para que de esta forma la diagonal sea la distancia entre el conductor y gnd
for k in range(cantidad_conductores):
    matriz_distancias[k][k]=coordenadasY[k]
#Introducir la altura de cada conductor en la diagonal de la matriz de distancias. Para que de esta forma la diagonal sea la distancia entre el conductor y gnd


#llamar a la función graficar
graficar2ConjuntosDeDatos(coordenadasX, gnd, "GND (tierra)", coordenadasY, "Conductor", "eje X (Metros [m])", "eje Y (Metros [m])", "Ubicacion en el plano cartesiano de los conductores ingresados","-","o")
#llamar a la función graficar
#---------------------------------------------------------------
#CALCULOS PREVIOS Y GRAFICA DEL SISTEMA DE CONDUCTORES INGRESADO
#---------------------------------------------------------------




#------------------------------
#CALCULOS DE PARAMETROS PROPIOS
#------------------------------
#Calculo de la resistencia de todos los conductores
matriz_resistencia=np.zeros(cantidad_conductores)  #crea una matriz de tamaño n, n siendo la cantidad de conductores
print(resistencia_general(radio_cable_guarda,resistividad_conductor,delta+radio_cable_guarda)[1])
for k in range(cantidad_conductores):
    if k>=cantidad_conductores-cantidad_cables_guarda:
        resistencia_k_cable=resistencia_general(radio_cable_guarda,resistividad_conductor,delta+radio_cable_guarda)[0]  #esto calculara la resistencia de los cables de guarda. Pero es muy importante entender que en el apartado de delta se coloca delta+radio_cable_guarda para que SIEMPRE sin importar la frecuencia de la linea de transmision la resistencia de los cables de guarda se calculen con la formula que no tienen en cuenta el efecto piel, porque estos no sufren este efecto
    else:
        resistencia_k_cable=resistencia_general(radio_conductores_fases,resistividad_conductor,delta)[0]  #esto calculara la resistencia de los cables de fase
    matriz_resistencia[k]=resistencia_k_cable

matriz_resistencia=1000*matriz_resistencia    #multiplicamos por 1000 para obtener [Ohms/km]
print("Resistencia de los conductores [Ohms/km]:")
print(matriz_resistencia)   
#Calculo de la resistencia de todos los conductores


#Calculo de la inductancia de propia y mutua de todos los conductores
matriz_inductancia=np.zeros((cantidad_conductores,cantidad_conductores),dtype=complex)  #crea una matriz de tamaño nxn, n siendo la cantidad de conductores
for k in range(cantidad_conductores-cantidad_cables_guarda):  #este bucle calcula la inductancia propia e introduce este valor en la diagonal de la matriz de inductancias de los conductores de fase
    matriz_inductancia[k][k]=inductancia_propia_aprox_efecto_piel_tierra_y_sin_efecto_piel(radio_conductores_fases,mu_medio,matriz_distancias[k][k],p,mu_conductor,resistividad_conductor,frecuencia)

for k in range(cantidad_cables_guarda):    #este bucle calcula la inductancia propia e introduce este valor en la diagonal de la matriz de inductancias de los cables de guarda
    matriz_inductancia[-cantidad_cables_guarda+k][-cantidad_cables_guarda+k]=inductancia_propia_aprox_efecto_piel_tierra_y_sin_efecto_piel(radio_cable_guarda,mu_medio,matriz_distancias[-cantidad_cables_guarda+k][-cantidad_cables_guarda+k],p,mu_conductor,resistividad_conductor,frecuencia)

for i in range(cantidad_conductores):  #Calculo de inductancias mutuas
    for j in range(cantidad_conductores):
        if i != j:
            coordenadas_para_la_imagen = coordenadas_conductores.copy() #hace que coordenadas_para_la_imagen tenga los mismos valores de la matriz que contiene las coordenadas d elos conductores
            coordenadas_para_la_imagen[j][1] = -coordenadas_para_la_imagen[j][1]  #reflejar respecto al suelo
            matriz_distancias_imagen = calcular_distancias(coordenadas_para_la_imagen, cantidad_conductores)   #calcula una matriz de distancias, pero ahora teniendo en cuenta las distancias entre la imagen y los otros conductores
            
            distancia_ij = matriz_distancias[i][j]                      # distancia entre conductor i y j
            distancia_imagen_ij = matriz_distancias_imagen[i][j]        # distancia entre conductor i y la imagen de j
            
            L_mutua = inductancia_mutua_aproximada(mu_medio, distancia_ij, distancia_imagen_ij)
            matriz_inductancia[i][j] = L_mutua

matriz_inductancia=1000*matriz_inductancia
#print("Matriz de inductancias [H/km]:")
#print(matriz_inductancia_metro)
#Calculo de la inductancia propia y mutua de todos los conductores


#Calculo de la matriz de capacitancia aprovechando la formula que la relaciona con la inductancia
matriz_capacitancia_metro=mu_medio*epsilon_medio*(np.linalg.inv(matriz_inductancia))
print("Matriz capacitancias [F/km]:")
print(matriz_capacitancia_metro)
#Calculo de la matriz de capacitancia aprovechando la formula que la relaciona con la inductancia
#------------------------------
#CALCULOS DE PARAMETROS PROPIOS
#------------------------------




#-----------------------------------------------------
#CREACION DE LA MATRIZ DE IMPEDANCIAS O MATRIZ DE FASE
#-----------------------------------------------------
matriz_fase=np.zeros((cantidad_conductores,cantidad_conductores),dtype=complex)  #crea una matriz de tamaño nxn, n siendo la cantidad de conductores

for k in range(cantidad_conductores):   #estos bucles van a recorrer cada elemento del array y van a calcular la impedancia de este
    for j in range(cantidad_conductores):
        if k==j:   #ya que matriz_resistencia_metro es un vector 1xn y matriz_inductancia_metro es una matriz nxn, no se pueden sumar. Pero hay que tener en cuenta que matriz_resistencia_metro solo se tiene en cuenta en la diagonal de la matriz
            impedancia=matriz_resistencia[k]+(1j*w*matriz_inductancia[k][j])
        else:
            impedancia=1j*w*matriz_inductancia[k][j]
        matriz_fase[k][j]=impedancia    #almacena el valor d ela impedancia en su ubicacion correspondiente en la matriz de fase
print("Matriz de fase [Ohms/km]:")
print(matriz_fase)
#-----------------------------------------------------
#CREACION DE LA MATRIZ DE IMPEDANCIAS O MATRIZ DE FASE
#-----------------------------------------------------




#--------------------------------------------------------------------------------------
#CREACION DE LA MATRIZ DE EQUIVALENTE PARA "ELIMINAR" EL EFECTO DE LOS CABLES DE GUARDA
#--------------------------------------------------------------------------------------
filas_y_columnas_para_matriz_principal = cantidad_conductores - cantidad_cables_guarda   #Variable que define el tamaño de la matriz de fases (sin guardas)
n_fase = filas_y_columnas_para_matriz_principal #Usamos n_fase para que el slicing sea más claro


#Sub-matriz principal
submatriz_principal = matriz_fase[0:n_fase, 0:n_fase]
#print(f"Sub-matriz principal {filas_y_columnas_para_matriz_principal}x{filas_y_columnas_para_matriz_principal} [Ohms/m]:")
#print(submatriz_principal)
#Sub-matriz principal


#Sub-matriz inferior
submatriz_inferior = matriz_fase[n_fase:, 0:n_fase]
#print(f"Sub-matriz inferior {cantidad_cables_guarda}x{filas_y_columnas_para_matriz_principal} [Ohms/m]:")
#print(submatriz_inferior)
#Sub-matriz inferior


#Sub-matriz lateral
submatriz_lateral = matriz_fase[0:n_fase, n_fase:]
#print(f"Sub-matriz lateral {filas_y_columnas_para_matriz_principal}x{cantidad_cables_guarda} [Ohms/m]:")
#print(submatriz_lateral)
#Sub-matriz lateral


#Sub-matriz cables de guarda (Z_gg)
submatriz_cables_guarda = matriz_fase[n_fase:, n_fase:]
#print(f"Sub-matriz cables de guarda {cantidad_cables_guarda}x{cantidad_cables_guarda} [Ohms/m]:")
#print(submatriz_cables_guarda)
#Sub-matriz cables de guarda


#Calculo de la matiz equivalente (Esta parte de tu código ya estaba bien)
matrices_para_calcular_matriz_equivalente = [submatriz_lateral, np.linalg.inv(submatriz_cables_guarda), submatriz_inferior]
matriz_equivalente = submatriz_principal - np.linalg.multi_dot(matrices_para_calcular_matriz_equivalente)
print(f"Matriz equivalente {filas_y_columnas_para_matriz_principal}x{filas_y_columnas_para_matriz_principal} [Ohms/km]")
print(matriz_equivalente)
#Calculo de la matiz equivalente
#--------------------------------------------------------------------------------------
#CREACION DE LA MATRIZ DE EQUIVALENTE PARA "ELIMINAR" EL EFECTO DE LOS CABLES DE GUARDA
#--------------------------------------------------------------------------------------




#----------------------------------
#CREACION DE LA MATRIZ DE SECUENCIA
#----------------------------------
#Generacion de la matriz de fortescue
T, T_inv=generar_matriz_fortescue(cantidad_circuitos_trifasicos)
#Generacion de la matriz de fortescue


#Calculo de la matriz de secuencia
matrices_para_calcular_matriz_secuencia=[T_inv,matriz_equivalente,T]
matriz_secuencia=np.linalg.multi_dot(matrices_para_calcular_matriz_secuencia)
print("Matriz de secuencia [Ohms/km]:")
print(matriz_secuencia)
#Calculo de la matriz de secuencia
#----------------------------------
#CREACION DE LA MATRIZ DE SECUENCIA
#----------------------------------