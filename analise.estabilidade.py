import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
from numpy.linalg import eig
import matplotlib.pyplot as plt

def modelo_estradiol(y, t, params):
    ( ke_1, ke_2, a_ed,d_o, k_12, k_21, d_1, d_2, k_o, k_1, d_c1, k_2, k_3, d_c2, K_lm, k_lb, d_b, p_cs, q_cd1, q_cd2, p_bs, q_bd, aed, kmax, Mmax, k_01, a_01, k_02, a_02, a_12, a_22, kpm, apm, apm1, dm, amb1, 
    kpb, apb, ae2, E2max, de2
    ) = params
    
    D, Mo, M1, M2, C1, C2, Cm, Cb, Mc, Mb, E2 = y
    
    # Phagocytosis rate

    def RD(t,D,a_ed):
      return (D/(a_ed+D))

    # Migration rate of M_o to the injury site

    def RM(t,Mo,M1,M2,kmax,Mmax):
      M = Mo+M1+M2
      return kmax*(1-M/Mmax)

    # Rate of differentiation of M_o into M_1

    def G1(t,k_01,a_01,C1):
      return k_01 * (C1/(a_01+C1))

    # Rate of differentiation of M_1 into M_0

    def G2(t,k_02,a_02,C2):
      return k_02 * (C2/(a_02+C2))

    # Inhibitory effect of the pro-inflammatory cytokine

    def H1(t,a_12,C2):
      return (a_12/(a_12+C2))

    # Inhibitory effect of the anti-inflammatory cytokine

    def H2(t,a_22,C2):
      return (a_22/(a_22+C2))

    # Inhibitory effect of the anti-inflammatory cytokine

    def Am(t,kpm,apm,apm1,C1):
      return kpm * ((apm*apm) + (apm1*C1)) / ((apm*apm) + (C1*C1))

    # Proliferation of C_m inhibited by c_1

    def F1(t,dm,amb1,C1):
      return dm * (amb1/(amb1+C1))

    # Proliferation of C_B inhibited by c_1

    def Ab(t,kpb,apb,C1):
      return kpb * (apb/(apb+C1))



    dDdt  = -RD(t,D,aed)*(ke_1*M1 + ke_2*M2)
    
    dModt = RM(t,Mo,M1,M2,kmax,Mmax) - G1(t,k_01,a_01,C1)*Mo - G2(t,k_02,a_02,C2)*Mo - (d_o*Mo)
    
    dM1dt  = G1(t,k_01,a_01,C1)*Mo + (k_21*M2) - (k_12*M1) - (d_1*M1)
    
    dM2dt  = G2(t,k_02,a_02,C2)*Mo 
    + (k_12*M1) - (k_21*M2) - (d_2*M2)
    dC1dt  = H1(t,a_12,C2)*(k_o*D + k_1*M1) - (d_c1*C1)*E2
    
    dC2dt  = H2(t,a_22,C2)*(k_2*M2 + k_3*Cm)*E2 - (d_c2*C2)
    
    dCmdt  = Am(t,kpm,apm,apm1,C1)*Cm*(1 - (Cm/K_lm)) - (F1(t,dm,amb1,C1)*Cm)
    
    dCbdt  = Ab(t,kpb,apb,C1)*Cb*(1 - (Cb/k_lb)) + (F1(t,dm,amb1,C1)*Cm) - (d_b*Cb)
    
    dMcdt  = (p_cs - q_cd1*Mc)*Cm - (q_cd2*Mc*Cb)
    
    dMbdt  = (p_bs - (q_bd*Mb))*Cb
    
    dE2dt  = ((E2max-E2)/E2max)*ae2*E2 - de2*E2

    return [
        dDdt, dModt, dM1dt, dM2dt, dC1dt, dC2dt, dCmdt, dCbdt, dMcdt, dMbdt, dE2dt
    ]

params = [
    3, 3, 4.71*10**6, 0.156, 0.075, 
    0.005, 0.121, 0.163, 0.0000005, 0.0000083, 0.2, 0.00000372, 0.0000007, 0.25,
    1000000.0, 1000000.0, 0.15,
    0.000003, 0.000003, 0.000002, 0.00000005, 0.00000005,
    4.71*10**6, 0.015, 6.0*10**5, 0.55, 0.01,
    0.3, 0.005, 0.025,
    0.1, 0.5, 3.162, 13.0, 1.0,
    0.1, 0.2202, 10.0,
    0.5, 0.019, 0.03
]

# === Condições iniciais (pós-menopausa como exemplo) ===
y0 = [500000, 4000, 0, 0, 100, 1, 1000, 0, 0, 0, 0.019]

#Ponto de equilíbrio

def steady_state(y_guess):
    return modelo_estradiol(y_guess, 0, params)

y_eq = fsolve(steady_state, y0)

print("Ponto de equilíbrio encontrado:")
print(y_eq)

#  Jacobiana 

def calc_jacobian(f, y_eq, params, eps=1e-8):
    n = len(y_eq)
    J = np.zeros((n, n))
    f0 = np.array(f(y_eq, 0, params))
    for i in range(n):
        y_perturb = y_eq.copy()
        y_perturb[i] += eps
        fi = np.array(f(y_perturb, 0, params))
        J[:, i] = (fi - f0) / eps
    return J

J = calc_jacobian(modelo_estradiol, y_eq, params)

print("\nMatriz Jacobiana no ponto de equilíbrio:")
print(J)

# autovalores

autovalores = eig(J)[0]

print("\nAutovalores da Jacobiana:")
print(autovalores)


if np.all(np.real(autovalores) < 0):
    print("\n O ponto de equilíbrio é estável, autovalores com parte real negativa.")
else:
    print("\n O ponto de equilíbrio é instável, algum autovalor com parte real positiva.")



    
