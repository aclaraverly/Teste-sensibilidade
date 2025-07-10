
# Import necessary libraries.
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint,solve_ivp
import post_processing as pp


# Note: eval is part of the Python standard library.
# More examples at the link: https://www.programiz.com/python-programming/methods/built-in/eval

"""Now, let's define the function that contains the equations and their respective parameters."""

# Define the function that contains the equations of the system to be evaluated.
def bonerepair(y,t, E2_max):
  # Parameters of the system of ten equations described above:
  ke_1 = 3.0
  ke_2 = 3.0
  a_ed = 4.71*10**6
  d_o = 0.156
  k_12 = 0.075
  k_21 = 0.005
  d_1 = 0.121
  d_2 = 0.163
  k_o = 0.0000005
  k_1 = 0.0000083
  d_c1 = 0.2 
  k_2 = 0.00000372
  k_3 = 0.0000007
  d_c2 = 0.25
  K_lm = 1000000.0
  k_lb = 1000000.0
  d_b = 0.15
  p_cs = 0.000003
  q_cd1 = 0.000003
  q_cd2 = 0.000002
  p_bs = 0.00000005
  q_bd = 0.00000005
  aed = 4.71*10**6
  kmax = 0.015
  Mmax = 6.0*10**5
  k_01 = 0.55
  a_01 = 0.01
  k_02 = 0.3
  a_02 = 0.005
  a_12 = 0.025
  a_22 = 0.1
  kpm = 0.5
  apm = 3.162
  apm1 = 13.0
  dm = 1.0
  amb1 = 0.1
  kpb = 0.2202
  apb = 10.0
  ae2 = 0.5 
  E2max = E2_max#0.019 
  de2 = 0.03  

  # Decouple to simplify the writing of the equations.
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
  
  # define the equations:
  dDdt  = -RD(t,D,aed)*(ke_1*M1 + ke_2*M2)
  dModt = RM(t,Mo,M1,M2,kmax,Mmax) - G1(t,k_01,a_01,C1)*Mo - G2(t,k_02,a_02,C2)*Mo - (d_o*Mo)
  dM1dt  = G1(t,k_01,a_01,C1)*Mo + (k_21*M2) - (k_12*M1) - (d_1*M1)
  dM2dt  = G2(t,k_02,a_02,C2)*Mo + (k_12*M1) - (k_21*M2) - (d_2*M2)
  dC1dt  = H1(t,a_12,C2)*(k_o*D + k_1*M1) - (d_c1*C1)*E2
  dC2dt  = H2(t,a_22,C2)*(k_2*M2 + k_3*Cm)*E2 - (d_c2*C2)
  dCmdt  = Am(t,kpm,apm,apm1,C1)*Cm*(1 - (Cm/K_lm)) - (F1(t,dm,amb1,C1)*Cm)
  dCbdt  = Ab(t,kpb,apb,C1)*Cb*(1 - (Cb/k_lb)) + (F1(t,dm,amb1,C1)*Cm) - (d_b*Cb)
  dMcdt  = (p_cs - q_cd1*Mc)*Cm - (q_cd2*Mc*Cb)
  dMbdt  = (p_bs - (q_bd*Mb))*Cb
  dE2dt  = ((E2max-E2)/E2max)*ae2*E2 - de2*E2

  # Recombines to return
  dydt = [dDdt, dModt, dM1dt, dM2dt, dC1dt, dC2dt, dCmdt, dCbdt, dMcdt, dMbdt, dE2dt]
  return dydt

def main():
  
  days = 100
  t = np.linspace(0, days, days)
  print(t)

  # condicoes iniciais
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
  print('Roda com o valor maximo...')
  # E2_max no valor maximo primeiro
  ys = odeint(bonerepair, yinit, t, args=(E2,))

  #D, Mo, M1, M2, C1, C2, Cm, Cb, Mc, Mb, E2
  # salva figuras
  #pp.plots(t,days,ys)

  #Atualiza condicao inicial do E2 para rodar a segunda vez
  E2 = 0.019
  yinit = np.array([D, Mo, M1, M2, C1, C2, Cm, Cb, Mc, Mb, E2])
  print('Roda com o valor minimo...')
  # E2_max com o valor minimo na segunda chamada
  ys_min = odeint(bonerepair, yinit, t, args=(E2,))
  pp.plots(t, days, ys)
  
  pp.plots_estradiol(t,days,ys,ys_min)

if __name__ == "__main__":
    main()