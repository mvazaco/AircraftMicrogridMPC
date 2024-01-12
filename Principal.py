import pyomo.environ as pyo
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Modelos import build_model_1, build_model_3
from Funciones import prepare_loads_data, extract_horizon_solution


cont = 0
horizon = 150

(P_L_crit_dic, P_L_HP_dic, P_L_LP_dic) = prepare_loads_data(cont, horizon)

mpc1 = build_model_1(horizon, P_L_crit_dic, P_L_HP_dic, P_L_LP_dic)

(t, P_in_aux, S_L_HP_aux, S_L_LP_aux, SOC_aux, P_ch_aux, P_disch_aux) = extract_horizon_solution(mpc1)
    
x = np.arange(0, 150)

P_L_crit = list(P_L_crit_dic.values())
P_L_HP = list(P_L_HP_dic.values())
P_L_LP = list(P_L_LP_dic.values())


plt.figure(figsize=(10, 5)) # POWER BALANCE

P_generated = []
for i in range(len(t)): P_generated.append( P_in_aux[i] + P_disch_aux[i] )
plt.subplot(211) 
plt.step(x, P_generated, color='green')
plt.xticks(np.arange(0, 160, 25))
plt.ylabel('Potencia generada [kW]')
plt.title('Balance de portencia', fontsize=20)
plt.grid(True)


P_consumed = []
for i in range(len(t)): P_consumed.append((-1)*(P_ch_aux[i] + np.multiply(S_L_HP_aux, P_L_HP)[i] + np.multiply(S_L_LP_aux, P_L_LP)[i] + P_L_crit[i]))
plt.subplot(212) 
plt.step(x, P_consumed, color='red')
plt.xlabel('Tiempo [min]')
plt.xticks(np.arange(0, 160, 25))
plt.ylabel('Potencia consumida [kW]')
plt.grid(True)

plt.show() ##################


plt.figure(figsize=(10, 7.5)) # POWER DEMAND / SHEDDING

plt.subplot(311)
plt.step(x, P_L_crit, label='Cargas críticas')
plt.xticks(np.arange(0, 160, 25))
plt.xlim(0)
plt.ylabel('Potencia[kW]')
plt.legend()
plt.grid(True)

plt.subplot(312)
plt.step(x, P_L_HP, label='Cargas de alta prioridad')
plt.step(x, np.multiply(P_L_HP, S_L_HP_aux), label='Con desconexión')
plt.xticks(np.arange(0, 160, 25))
plt.xlim(0)
plt.ylabel('Potencia[kW]')
plt.legend(loc = 'upper left')
plt.grid(True)

plt.subplot(313)
plt.step(x, P_L_LP, label='Cargas de baja prioridad')
plt.step(x, np.multiply(P_L_LP, S_L_LP_aux), label='Con desconexión')
plt.xlabel('Tiempo [min]')
plt.xticks(np.arange(0, 160, 25))
plt.xlim(0)
plt.ylabel('Potencia[kW]')
plt.legend()
plt.grid(True)

plt.show() ######################

plt.figure(figsize=(10,5)) # SOC
plt.plot(x, SOC_aux)
plt.xlabel('Tiempo [min]')
plt.xticks(np.arange(0, 160, 25))
plt.xlim(0,150)
plt.ylim(0.2, 1)
plt.grid(True)
plt.title('SOC', fontsize=20)

LO = []
for i in range(len(x)): LO.append(0.3)
HI = []
for i in range(len(x)): HI.append(0.9)
plt.plot(x, HI, color='red')
plt.plot(x, LO, color='red')

plt.show() ######################

plt.figure(figsize=(10,5)) # CHARGE/DISCHARGE

plt.step(x, P_ch_aux, color='green', label='Potencia de carga')
plt.xticks(np.arange(0, 160, 25))
plt.xlim(0,150)
plt.grid(True)


P_disch_aux_neg = np.multiply(P_disch_aux, -1)
plt.step(x, P_disch_aux_neg, color='red', label='Potencia de descarga')
plt.xticks(np.arange(0, 160, 25))
plt.xlim(0,150)
plt.legend(loc='upper right', fontsize='large')

plt.show() #######################