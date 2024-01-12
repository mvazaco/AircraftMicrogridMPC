import pandas as pd
import pyomo.environ as pyo
import numpy as np
import matplotlib.pyplot as plt

mpc = pyo.ConcreteModel()

############################# PARAMETERS #################################

mpc.N = pyo.Param(initialize=5)
mpc.T = pyo.Param(initialize=15)
mpc.T_s = pyo.Param(initialize=0.0166)

mpc.P_in_max  = pyo.Param(initialize=4)

mpc.B_cap = pyo.Param(initialize=4)
mpc.SOC_0 = pyo.Param(initialize=0.3)
mpc.LO = pyo.Param(initialize=0.3)
mpc.HI = pyo.Param(initialize=0.9)

mpc.P_ch_max    = pyo.Param(initialize=4)
mpc.P_disch_max = pyo.Param(initialize=4)
mpc.nu_ch    = pyo.Param(initialize=0.9)
mpc.nu_disch = pyo.Param(initialize=0.9)

mpc.N_L = pyo.Param(initialize=2)

mpc.omega_S       = pyo.Param(initialize=5)
mpc.omega_SOC     = pyo.Param(initialize=16)
mpc.omega_P_batt  = pyo.Param(initialize=2)
mpc.omega_P_ch    = pyo.Param(initialize=3.1)
mpc.omega_delta   = pyo.Param(initialize=3.47)
mpc.omega_Delta   = pyo.Param(initialize=1e+05)
mpc.omega_P_disch = pyo.Param(initialize=5.1)


mpc.k = pyo.Set(initialize = range(mpc.N()))
mpc.i = pyo.Set(initialize = range(mpc.N_L()))

mpc.gamma_L = pyo.Param(mpc.i, initialize=[2,1])
mpc.S_L_0   = pyo.Param(mpc.i, initialize=[1,0])
# mpc.S_L_0.pprint()
                                
critical_loads_data = [1, 1, 1, 4, 4, 3, 2, 2, 2, 2, 2, 3, 3, 1, 0];            
HP_loads_data   = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 3, 3, 2, 1];
LP_loads_data   = [2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3, 3, 1, 0];

# plt.plot(critical_loads_data)
# plt.show()

critical_loads_window = {}
for i,v in enumerate(critical_loads_data):
    if i < pyo.value(mpc.N()):
        critical_loads_window[i] = v
    
HP_loads_window = {}
for i,v in enumerate(HP_loads_data):
    if i < pyo.value(mpc.N()):
        HP_loads_window[i] = v
    
LP_loads_window = {}
for i,v in enumerate(LP_loads_data):
    if i < pyo.value(mpc.N()):
        LP_loads_window[i] = v

mpc.P_L_crit = pyo.Param(mpc.k, initialize = critical_loads_window)
mpc.P_L_HP = pyo.Param(mpc.k, initialize = HP_loads_window)
mpc.P_L_LP = pyo.Param(mpc.k, initialize = LP_loads_window)

############################# VARIABLES #################################

mpc.P_in    = pyo.Var(mpc.k, bounds=(0,4))

mpc.S_L     = pyo.Var(mpc.i,mpc.k, domain=pyo.Binary)
mpc.lineal  = pyo.Var(mpc.i,mpc.k, within=pyo.NonNegativeReals)

mpc.SOC         = pyo.Var(mpc.k, bounds=(mpc.LO,mpc.HI))

mpc.P_ch        = pyo.Var(mpc.k, bounds=(0,mpc.P_ch_max))
mpc.P_disch     = pyo.Var(mpc.k, bounds=(0,mpc.P_disch_max))

mpc.dseta_disch = pyo.Var(mpc.k, domain=pyo.Binary)
mpc.dseta_ch    = pyo.Var(mpc.k, domain=pyo.Binary)


############################# OBJECTIVE FUNCTION #################################

def obj1_expression(mpc):
    
    J_S_L = sum(mpc.gamma_L[i] * (1 - mpc.S_L[i,k]) / (mpc.N * (mpc.gamma_L[0] + mpc.gamma_L[1])) for k in mpc.k for i in mpc.i)
    
    J_delta_L = sum( mpc.lineal[i,k] / (mpc.N * mpc.N_L) for k in mpc.k for i in mpc.i)
    
    J_P_ch = sum((mpc.P_ch_max - mpc.P_ch[k]) / (mpc.N * mpc.P_ch_max) for k in mpc.k)
      
    J_P_disch = sum(mpc.P_disch[k] / (mpc.N * mpc.P_disch_max) for k in mpc.k)

    return mpc.omega_S*J_S_L + mpc.omega_delta*J_delta_L + mpc.omega_P_ch*J_P_ch + mpc.omega_P_disch*J_P_disch
    # return mpc.omega_S*J_S_L + mpc.omega_P_ch*J_P_ch + mpc.omega_P_disch*J_P_disch

mpc.obj1 = pyo.Objective(expr = obj1_expression)

############################# CONSTRAINTS #################################

def power_balance(mpc, k):
    return mpc.P_in[k] - mpc.P_ch[k] + mpc.P_disch[k] - mpc.S_L[0,k]*mpc.P_L_HP[k] - mpc.S_L[1,k]*mpc.P_L_LP[k] - mpc.P_L_crit[k] == 0
mpc.power_balance_constraint = pyo.Constraint(mpc.k, rule=power_balance)

def SOC_def(mpc, k):
    if(k==0):
        return mpc.SOC_0 + mpc.T_s*mpc.nu_ch*mpc.P_ch[k]/mpc.B_cap - mpc.T_s*mpc.P_disch[k]/(mpc.nu_disch*mpc.B_cap) == mpc.SOC[k]
    else:
        return mpc.SOC[k-1] + mpc.T_s*mpc.nu_ch*mpc.P_ch[k]/mpc.B_cap - mpc.T_s*mpc.P_disch[k]/(mpc.nu_disch*mpc.B_cap) == mpc.SOC[k]
mpc.SOC_constraint = pyo.Constraint(mpc.k, rule=SOC_def)

# def battery_mode(mpc, k):
#     return mpc.P_ch[k] * mpc.P_disch[k] <= 1e-03
# mpc.SOC_constraint = pyo.Constraint(mpc.k, rule=battery_mode)

def battery_mode(mpc, k):
    return  mpc.dseta_ch[k] + mpc.dseta_disch[k] == 1
mpc.battery_mode_constraint = pyo.Constraint(mpc.k, rule=battery_mode)

def charging_mode(mpc, k):
    return  mpc.P_ch[k] <= mpc.dseta_ch[k]*mpc.P_ch_max
mpc.charging_mode_constraint = pyo.Constraint(mpc.k, rule=charging_mode)

def discharging_mode(mpc, k):
    return  mpc.P_disch[k] <=  mpc.dseta_disch[k]*mpc.P_disch_max
mpc.discharging_mode_constraint = pyo.Constraint(mpc.k, rule=discharging_mode)

def lineal_pos (mpc, i, k):
    if(k==0):
        return mpc.S_L[i,k] - mpc.S_L_0[i] <= mpc.lineal[i,k]
    else:
        return mpc.S_L[i,k] - mpc.S_L[i,k-1] <= mpc.lineal[i,k]
mpc.lineal_pos_constraint = pyo.Constraint(mpc.i,mpc.k, rule=lineal_pos)
    
def lineal_neg (mpc, i, k):
    if(k==0):
        return mpc.S_L_0[i] - mpc.S_L[i,k] <= mpc.lineal[i,k]
    else:
        return mpc.S_L[i,k-1] - mpc.S_L[i,k]  <= mpc.lineal[i,k]
mpc.lineal_neg_constraint = pyo.Constraint(mpc.i,mpc.k, rule=lineal_neg)



opt = pyo.SolverFactory('cplex')
results = opt.solve(mpc)
print('Solver Status: ', results.solver.termination_condition)



# mpc.solutions.store_to(results)
# results.write(filename='results.json', format='json')

# mpc.display()
mpc.pprint()

# for i in range(0,3):
# for t in range(0,24):
# BuyGrid[i,t] = M1.BuyGrid[i,t].value
