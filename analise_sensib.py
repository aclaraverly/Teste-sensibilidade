
from bone_estradiol import bonerepair

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import SALib
from SALib.analyze import sobol
from SALib.sample import saltelli


problem = {
    'num_vars': 2,
    'names': ['ke_1', 'ke_2'],
    'bounds': [[2.5, 3.5], [2.5, 3.5]] 
}

param_values = saltelli.sample(problem, 1024)

Y = []

D  = 5.0*10**5
Mo = 4000.0
M1 = 0.0
M2 = 0.0
C1 = 100.0
C2 = 1.0
Cm = 1000.0 
Cb = 0.0 
Mc = 0.0
Mb = 0.0
E2 = 0.060


yinit = np.array([D, Mo, M1, M2, C1, C2, Cm, Cb, Mc, Mb, E2])

for params in param_values:
    ke_1, ke_2 = params
    sol = solve_ivp(bonerepair, [0, 10], yinit, args=(ke_1,ke_2,))
    Y.append(sol.y[7][-1])  #qual a variavel sera influenciada pelo parametro

Y = np.array(Y)  

Si = sobol.analyze(problem, Y, print_to_console=True)  # Método Sobol, que decompõe a variância da saída em contribuições de cada parâmetro.
#'print_to_console=True' exibe os resultados no console.

df_sobol = pd.DataFrame({
    'Parameter': problem['names'],
    'S1': Si['S1'],          # S1: Índice de sensibilidade de primeira ordem (contribuição direta de um parâmetro);
    'ST': Si['ST'],          # ST: Índice de sensibilidade total (considera também interações entre parâmetros);
    'S1_conf': Si['S1_conf'],  # Intervalos de confiança dos índices.
    'ST_conf': Si['ST_conf']   # Intervalos de confiança dos índices.
})


print(df_sobol) # Permite ver quais parâmetros mais influenciam a variável de saída D.


# Influência de cada parâmetro na saída, com barras de erro representando os intervalos de confiança.

plt.figure(figsize=(8, 6))
plt.bar(df_sobol['Parameter'], df_sobol['S1'], yerr=df_sobol['S1_conf'], capsize=5)
plt.xlabel('Parâmetros')
plt.ylabel('Índice de Sensibilidade S1')
plt.title('Análise de Sensibilidade - Sobol (S1)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()


plt.figure(figsize=(8, 6))
plt.bar(df_sobol['Parameter'], df_sobol['ST'], yerr=df_sobol['ST_conf'], capsize=5, color='orange')
plt.xlabel('Parâmetros')
plt.ylabel('Índice de Sensibilidade Total (ST)')
plt.title('Análise de Sensibilidade - Sobol (ST)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

df_sobol.to_csv('resultados_sobol.csv', index=False, encoding='utf-8') 