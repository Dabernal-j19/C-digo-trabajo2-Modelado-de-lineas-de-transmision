import numpy as np

np.set_printoptions(linewidth=400)  # Esto ajusta la terminal para mostrar toda la fila completa

# Aqui colocar el array resultado arrojado por Python
matriz_fase_python = np.array([
    [0.25101338+0.85016545j, 0.21151553j, 0.16566683j, 0.13509548j],
    [0.21151553j, 0.25101338+0.85904609j, 0.21990025j, 0.16673416j],
    [0.16566683j, 0.21990025j, 0.25101338+0.86699008j, 0.21651954j],
    [0.13509548j, 0.16673416j, 0.21651954j, 0.53398803+0.90399269j]
], dtype=complex)

# Aqui colocar el array resultado arrojado por ATP
matriz_fase_ATP = np.array([
    [7.374047e-4+6.648553e-1j, 1.152456e-6+2.2115171e-1j, 1.088649e-6+1.656682e-1j, 1.020609e-6+1.350968e-1j],
    [1.152456e-6+2.2115171e-1j, 7.372686e-4+6.737358e-1j, 1.031188e-6+2.199017e-1j, 9.699475e-7+1.667355e-1j],
    [1.088649e-6+1.656682e-1j, 1.031188e-6+2.199017e-1j, 7.371598e-4+6.816797e-1j, 9.243263e-7+2.165209e-1j],
    [1.020609e-6+1.350968e-1j, 9.699475e-7+1.667355e-1j, 9.243263e-7+2.165209e-1j, 7.370548e-4+7.186823e-1j]
], dtype=complex)

# Tamaño de la matriz
tamaño = np.shape(matriz_fase_python)
print("Tamaño de las matrices:", tamaño)

# Convertir cada elemento a forma polar (módulo)
matriz_python_modulo = np.abs(matriz_fase_python)
matriz_ATP_modulo = np.abs(matriz_fase_ATP)

# Crear matriz vacía para errores
matriz_errores = np.zeros(tamaño)

# Calcular errores en porcentaje
for i in range(tamaño[0]):
    for j in range(tamaño[1]):
        if matriz_python_modulo[i][j] != 0:
            error = (abs(matriz_python_modulo[i][j] - matriz_ATP_modulo[i][j]) / matriz_python_modulo[i][j]) * 100
        else:
            error = 0
        matriz_errores[i][j] = error

print("\nLos errores de cada impedancia de la matriz de fase en % son:")
print(np.round(matriz_errores, 4))
