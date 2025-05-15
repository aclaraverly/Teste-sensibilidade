import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from bone_estradiol import bonerepair

# Tempo de simulação
dias = 360
t = np.linspace(0, dias, 361) # 0 a 360 dias.

# Condições iniciais, vetor com as onze variáveis.
yinit = np.array([
    500000,  # D
    4000.0, # Mo
    0.0,    # M1
    0.0,    # M2
    100.0,  # C1
    1.0,    # C2
    1000.0, # Cm
    0.0,    # Cb <- alvo da análise
    0.0,    # Mc
    0.0,    # Mb
    0.060   # E2
])

# Parâmetros do modelo, argumentos para as EDO's
params_base = {
    'ke_1': 3.0, 'ke_2': 3.0, 'a_ed': 4.71e6, 'd_o': 0.156, 'k_12': 0.075, 'k_21': 0.005,
    'd_1': 0.121, 'd_2': 0.163, 'k_o': 5e-7, 'k_1': 8.3e-6, 'd_c1': 0.2, 'k_2': 3.72e-6,
    'k_3': 7e-7, 'd_c2': 0.25, 'K_lm': 1e6, 'k_lb': 1e6, 'd_b': 0.15, 'p_cs': 3e-6,
    'q_cd1': 3e-6, 'q_cd2': 2e-6, 'p_bs': 5e-8, 'q_bd': 5e-8, 'aed': 4.71e6, 'kmax': 0.015,
    'Mmax': 6e5, 'k_01': 0.55, 'a_01': 0.01, 'k_02': 0.3, 'a_02': 0.005, 'a_12': 0.025,
    'a_22': 0.1, 'kpm': 0.5, 'apm': 3.162, 'apm1': 13.0, 'dm': 1.0, 'amb1': 0.1,
    'kpb': 0.2202, 'apb': 10.0, 'ae2': 0.5, 'E2max': 0.019, 'de2': 0.03
}

values = np.linspace(0.1, 1.1, 10)
prm = []
# Armazenar resultados de Cb para cada parametro
Cb_all = []

# modelo é simulado 10 vezes, cada uma com um valor diferente de cada parametro, para ver como isso afeta a variável Cb.
i=0
for v in values:
    params = params_base.copy() # cria uma cópia evitando alterar o dicionário original a cada iteração.
    params['de2'] = params['de2']*v # substitui por um novo valor.
    prm.append(params['de2'])
    i=i+1
    p_vals = list(params.values()) # Extrai apenas os valores numéricos dos parâmetros e os converte em uma lista.
    sol = odeint(bonerepair, yinit, t, args=tuple(p_vals))  
    Cb_all.append(sol[:, 7])  # Cb 
    
#print(prm)

Cb_all = np.array(Cb_all)

#print(Cb_all)

# Estatísticas temporais
Cb_mean = np.mean(Cb_all, axis=0)  #media temporal
Cb_std = np.std(Cb_all, axis=0)    #desvio


plt.figure(figsize=(10, 6))
for i, Cb in enumerate(Cb_all):
    plt.plot(t, Cb, alpha=0.3, label=f'de2 = {prm[i]:.4f}')
plt.plot(t, Cb_mean, color='black', label='Média', linewidth=2)
plt.xlabel('Dias')
plt.ylabel('Cb (Osteoblastos)')
plt.title('Sensibilidade de Cb ao longo do tempo em relação a de2')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



      

