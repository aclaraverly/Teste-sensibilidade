
from bone_estradiol import bonerepair

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import SALib
from SALib.analyze import sobol
from SALib.sample import saltelli


problem = {
    'num_vars': 41,
    'names': ['ke_1', 'ke_2', 'a_ed', 'd_o', 'k_12', 'k_21', 'd_1', 'd_2', 'k_o', 'k_1', 'd_c1', 'k_2', 'k_3', 'd_c2', 'K_lm', 'k_lb', 'd_b', 'p_cs', 'q_cd1', 'q_cd2', 'p_bs', 'q_bd', 'aed', 'kmax', 'Mmax', 'k_01', 'a_01', 'k_02', 'a_02', 'a_12', 'a_22', 'kpm', 'apm', 'apm1', 'dm', 'amb1', 'kpb', 'apb', 'ae2', 'E2max', 'de2'],
    'bounds': [[1.5, 4.5], [1.5, 4.5], [2355000, 7065000], [0.078, 0.234], [0.0375, 0.1125], [0.0025, 0.0075], [0.0605, 0.1815], [0.0815, 0.2445], [0.00000025, 0.00000075], [0.00000415, 0.00001245], [0.1, 0.3], [0.00000186, 0.00000558], [0.00000035, 0.00000105], [0.125, 0.375], [500000, 1500000], [500000, 1500000], [0.075, 0.225], [0.0000015, 0.0000045], [0.0000015, 0.0000045], [0.000001, 0.000003], [0.000000025, 0.000000075], [0.000000025, 0.000000075], [2355000, 7065000], [0.0075, 0.0225], [300000, 900000], [0.275, 0.825], [0.005, 0.015], [0.15, 0.45], [0.0025, 0.0075], [0.0125, 0.0375], [0.05, 0.15], [0.25, 0.75], [1.581, 4.743], [6.5, 19.5], [0.5, 1.5], [	0.05, 0.15], [0.1101, 0.3303], [5.0, 15.0], [0.25, 0.75], [0.0095, 0.0285], [0.015, 0.045]] 
}  # 50% -/+

param_values = saltelli.sample(problem, 2048)

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
    ke_1, ke_2, a_ed, d_o, k_12, k_21, d_1, d_2, k_o, k_1, d_c1, k_2, k_3, d_c2, K_lm, k_lb, d_b, p_cs, q_cd1, q_cd2, p_bs, q_bd, aed, kmax, Mmax, k_01, a_01, k_02, a_02, a_12, a_22, kpm, apm, apm1, dm, amb1, kpb, apb, ae2, E2max, de2 = params
    sol = solve_ivp(bonerepair, [0, 10], yinit, args=(ke_1, ke_2, a_ed, d_o, k_12, k_21, d_1, d_2, k_o, k_1, d_c1, k_2, k_3, d_c2, K_lm, k_lb, d_b, p_cs, q_cd1, q_cd2, p_bs, q_bd, aed, kmax, Mmax, k_01, a_01, k_02, a_02, a_12, a_22, kpm, apm, apm1, dm, amb1, kpb, apb, ae2, E2max, de2,))
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

plt.figure(figsize=(24, 18))
plt.bar(df_sobol['Parameter'], df_sobol['S1'], yerr=df_sobol['S1_conf'], capsize=5)
plt.xlabel('Parâmetros')
plt.ylabel('Índice de Sensibilidade S1')
plt.title('Análise de Sensibilidade - Sobol (S1)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()
plt.savefig()


plt.figure(figsize=(24, 18))
plt.bar(df_sobol['Parameter'], df_sobol['ST'], yerr=df_sobol['ST_conf'], capsize=5, color='orange')
plt.xlabel('Parâmetros')
plt.ylabel('Índice de Sensibilidade Total (ST)')
plt.title('Análise de Sensibilidade - Sobol (ST)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()
plt.savefig()

df_sobol.to_csv('resultados_sobol.csv', index=False, encoding='utf-8') 
