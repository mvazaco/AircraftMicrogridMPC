import pyomo.environ as pyo
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

# def plot_power_balance(x, horizon, P_L_HP, P_L_LP, P_L_crit, P_in, S_L_HP, S_L_LP, SOC, P_ch, P_disch):
    
#     plt.figure(figsize=(10, 5)) # POWER BALANCE

#     P_generated = []
#     for i in range(150): P_generated.append( P_in[i] + P_disch[i] )
#     plt.subplot(211) 
#     plt.step(x, P_generated, color='green')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.ylabel('Potencia generada [kW]')
#     plt.title('Balance de portencia (horizon = %i)' %horizon , fontsize=20)
#     plt.grid(True)



#     P_consumed = []
#     for i in range(150): P_consumed.append((-1)*(P_ch[i] + np.multiply(S_L_HP, P_L_HP)[i] + np.multiply(S_L_LP, P_L_LP)[i] + P_L_crit[i]))
#     plt.subplot(212) 
#     plt.step(x, P_consumed, color='red')
#     plt.xlabel('Tiempo [min]')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.ylabel('Potencia consumida [kW]')
#     plt.grid(True)

#     plt.show()
    
#     return

def plot_power_balance(x, horizon, P_L_HP, P_L_LP, P_L_crit, P_in, S_L, P_ch, P_disch):
    
    plt.figure(figsize=(10, 5)) # POWER BALANCE

    P_generated = []
    for i in range(150): P_generated.append( P_in[i] + P_disch[i] )
    plt.subplot(211) 
    plt.step(x, P_generated, color='green')
    plt.xticks(np.arange(0, 160, 25))
    plt.ylabel('Potencia generada [kW]')
    plt.title('Balance de portencia (horizonte = %i)' %horizon , fontsize=20)
    plt.grid(True)



    P_consumed = []
    for i in range(150): P_consumed.append((-1)*(P_ch[i] + np.multiply(S_L[0], P_L_HP)[i] + np.multiply(S_L[1], P_L_LP)[i] + P_L_crit[i]))
    plt.subplot(212) 
    plt.step(x, P_consumed, color='red')
    plt.xlabel('Tiempo [min]')
    plt.xticks(np.arange(0, 160, 25))
    plt.ylabel('Potencia consumida [kW]')
    plt.grid(True)

    plt.show()
    
    return
    

# def plot_demand_shedding(x, P_L_HP, P_L_LP, P_L_crit, P_in, S_L_HP, S_L_LP, SOC, P_ch, P_disch):

#     plt.figure(figsize=(10, 7.5)) # POWER DEMAND / SHEDDING

#     plt.subplot(311)
#     plt.step(x, P_L_crit, label='Cargas críticas')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.xlim(0)
#     plt.ylabel('Potencia[kW]')
#     plt.legend()
#     plt.grid(True)

#     plt.subplot(312)
#     plt.step(x, P_L_HP, label='Cargas de alta prioridad')
#     plt.step(x, np.multiply(P_L_HP, S_L_HP), label='Con desconexión')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.xlim(0)
#     plt.ylabel('Potencia[kW]')
#     plt.legend(loc = 'upper left')
#     plt.grid(True)

#     plt.subplot(313)
#     plt.step(x, P_L_LP, label='Cargas de baja prioridad')
#     plt.step(x, np.multiply(P_L_LP, S_L_LP), label='Con desconexión')
#     plt.xlabel('Tiempo [min]')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.xlim(0)
#     plt.ylabel('Potencia[kW]')
#     plt.legend()
#     plt.grid(True)

#     plt.show()
    
#     return

def plot_demand_shedding(x, horizon, P_L_HP, P_L_LP, P_L_crit, P_in, S_L, P_ch, P_disch):

    plt.figure(figsize=(10, 7.5)) # POWER DEMAND / SHEDDING

    plt.subplot(311)
    plt.step(x, P_L_crit, label='Cargas críticas')
    plt.xticks(np.arange(0, 160, 25))
    plt.xlim(0)
    plt.ylabel('Potencia[kW]')
    plt.legend()
    plt.grid(True)
    plt.title('horizonte = %i' %horizon , fontsize=20)

    plt.subplot(312)
    plt.step(x, P_L_HP, label='Cargas de alta prioridad')
    plt.step(x, np.multiply(P_L_HP, S_L[0]), label='Con desconexión')
    plt.xticks(np.arange(0, 160, 25))
    plt.xlim(0)
    plt.ylabel('Potencia[kW]')
    plt.legend(loc = 'upper left')
    plt.grid(True)

    plt.subplot(313)
    plt.step(x, P_L_LP, label='Cargas de baja prioridad')
    plt.step(x, np.multiply(P_L_LP, S_L[1]), label='Con desconexión')
    plt.xlabel('Tiempo [min]')
    plt.xticks(np.arange(0, 160, 25))
    plt.xlim(0)
    plt.ylabel('Potencia[kW]')
    plt.legend()
    plt.grid(True)

    plt.show()
    
    return


# def plot_SOC(x, P_L_HP, P_L_LP, P_L_crit, P_in, S_L_HP, S_L_LP, SOC, P_ch, P_disch):
    
#     plt.figure(figsize=(10,5)) # SOC
#     plt.plot(x, SOC)
#     plt.xlabel('Tiempo [min]')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.xlim(0,150)
#     plt.ylim(0.2, 1)
#     plt.grid(True)
#     plt.title('SOC', fontsize=20)

#     LO = []
#     for i in range(len(x)): LO.append(0.3)
#     HI = []
#     for i in range(len(x)): HI.append(0.9)
#     plt.plot(x, HI, color='red')
#     plt.plot(x, LO, color='red')

#     plt.show()
    
#     return

def plot_SOC(x, SOC):
    
    plt.figure(figsize=(10,5)) # SOC
    plt.plot(x, SOC)
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

    plt.show()
    
    return


# def plot_charge_discharge(x, P_L_HP, P_L_LP, P_L_crit, P_in, S_L_HP, S_L_LP, SOC, P_ch, P_disch):

#     plt.figure(figsize=(10,5)) # CHARGE/DISCHARGE

#     plt.step(x, P_ch, color='green', label='Potencia de carga')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.xlim(0,150)
#     plt.grid(True)


#     P_disch_neg = np.multiply(P_disch, -1)
#     plt.step(x, P_disch_neg, color='red', label='Potencia de descarga')
#     plt.xticks(np.arange(0, 160, 25))
#     plt.xlim(0,150)
#     plt.legend(loc='upper right', fontsize='large')

#     plt.show()
    
#     return

def plot_charge_discharge(x, horizon, P_ch, P_disch):

    plt.figure(figsize=(10,5)) # CHARGE/DISCHARGE

    plt.step(x, P_ch, color='green', label='Potencia de carga')
    plt.xticks(np.arange(0, 160, 25))
    plt.xlim(0,150)
    plt.grid(True)


    P_disch_neg = np.multiply(P_disch, -1)
    plt.step(x, P_disch_neg, color='red', label='Potencia de descarga')
    plt.xticks(np.arange(0, 160, 25))
    plt.xlim(0,150)
    plt.legend(loc='upper right', fontsize='large')
    plt.title('horizonte = %i' %horizon , fontsize=20)

    plt.show()
    
    return


def plot_G_S_L_h(horizons, G_S_L_h):

    plt.figure(figsize=(10,5))

    plt.plot(horizons, G_S_L_h, label='Obj1')
    plt.xlabel('Horizonte')
    plt.ylabel('G_S_L_h')
    plt.grid(True)
    plt.legend(fontsize='large')
    
    plt.scatter(horizons, G_S_L_h, marker='.')
    
    plt.show()
    
    return


















