import pandas as pd
import pyomo.environ as pyo
import numpy as np
import matplotlib.pyplot as plt


################################################ PYOMO ############################################

def build_model(horizon, P_L_crit, P_L_HP, P_L_LP):
    
    mpc = pyo.ConcreteModel()
    
    ####################### SACALAR PARAMETERS #####################
    
    mpc.N = pyo.Param(initialize=horizon)
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
    mpc.omega_delta   = pyo.Param(initialize=34.7)
    mpc.omega_Delta   = pyo.Param(initialize=1e+05)
    mpc.omega_P_disch = pyo.Param(initialize=5.1)
    
    ############################ SETS ##############################
    
    mpc.k = pyo.Set(initialize = range(mpc.N()))
    mpc.i = pyo.Set(initialize = range(mpc.N_L()))
    
    #################### INDEXED PARAMETERS ########################
    
    mpc.gamma_L = pyo.Param(mpc.i, initialize=[2,1])
    mpc.S_L_0   = pyo.Param(mpc.i, initialize=[1,0])
                                    
    ######################### VARIABLES ############################
    
    mpc.P_in    = pyo.Var(mpc.k, bounds=(0,4))
    
    mpc.S_L     = pyo.Var(mpc.i,mpc.k, domain=pyo.Binary)
    mpc.lineal  = pyo.Var(mpc.i,mpc.k, within=pyo.NonNegativeReals)
    
    mpc.SOC     = pyo.Var(mpc.k, bounds=(mpc.LO,mpc.HI))
    
    mpc.P_ch    = pyo.Var(mpc.k, bounds=(0,4))
    mpc.P_disch = pyo.Var(mpc.k, bounds=(0,4))
    
    
    ###################### OBJECTIVE FUNCTION ######################
    
    def obj1_expression(mpc):
        
        J_S_L = sum(mpc.gamma_L[i] * (1 - mpc.S_L[i,k]) / (mpc.N * (mpc.gamma_L[0]+mpc.gamma_L[1])) for k in mpc.k for i in mpc.i)
        
        J_delta_L = sum( mpc.lineal[i,k] / (mpc.N * mpc.N_L) for k in mpc.k for i in mpc.i)

                
        J_P_ch = sum((mpc.P_ch_max() - mpc.P_ch[k]) / (mpc.N() * mpc.P_ch_max()) for k in mpc.k)
          
        J_P_disch = sum(mpc.P_disch[k] / (mpc.N * mpc.P_disch_max) for k in mpc.k)
    
        return mpc.omega_S*J_S_L + mpc.omega_delta*J_delta_L + mpc.omega_P_ch*J_P_ch + mpc.omega_P_disch*J_P_disch
        # return mpc.omega_S*J_S_L + mpc.omega_P_ch*J_P_ch + mpc.omega_P_disch*J_P_disch
    
    mpc.obj1 = pyo.Objective(expr = obj1_expression)
    
    ############################# CONSTRAINTS ######################
    
    def power_balance(mpc, k):
        return mpc.P_in[k] - mpc.P_ch[k] + mpc.P_disch[k] - mpc.S_L[0,k]*P_L_HP[k] - mpc.S_L[1,k]*P_L_LP[k] - P_L_crit[k] == 0
    mpc.power_balance_constraint = pyo.Constraint(mpc.k, rule=power_balance)
    
    def SOC(mpc, k):
        if(k==0):
            return mpc.SOC_0 + mpc.T_s*mpc.nu_ch*mpc.P_ch[k]/mpc.B_cap - mpc.T_s*mpc.P_disch[k]/(mpc.nu_disch*mpc.B_cap) == mpc.SOC[k]
        else:
            return mpc.SOC[k-1] + mpc.T_s*mpc.nu_ch*mpc.P_ch[k]/mpc.B_cap - mpc.T_s*mpc.P_disch[k]/(mpc.nu_disch*mpc.B_cap) == mpc.SOC[k]
    mpc.SOC_constraint = pyo.Constraint(mpc.k, rule=SOC)
    
    def battery_mode(mpc, k):
        return mpc.P_ch[k] * mpc.P_disch[k] < 0
    
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
    
    ############################# SOLVER ###########################
    
    
    opt = pyo.SolverFactory('cplex')
    results = opt.solve(mpc)
    
    mpc.pprint()
    print('Solver Status: ', results.solver.termination_condition)

    # mpc.solutions.store_to(results)
    # results.write(filename='results.json', format='json')
    # mpc.display()
    
    return mpc
    

######################################### LOADS DATA ##################################################

# READ CSV LOADS DEMAND
critical_loads_dic = pd.read_csv ('D:\\Critical_loads.csv', names = ['critical loads'])
HP_loads_dic       = pd.read_csv ('D:\\HP_loads.csv', names = ['HP loads'])
LP_loads_dic       = pd.read_csv ('D:\\LP_loads.csv', names = ['LP loads'])

#PLOT LOADS DEMAND
critical_loads_np_array = critical_loads_dic["critical loads"][0:151].to_numpy()
HP_loads_np_array       = HP_loads_dic["HP loads"][0:151].to_numpy()
LP_loads_np_array       = LP_loads_dic["LP loads"][0:151].to_numpy()

load_profile, axs = plt.subplots(3, 1, figsize=(5,5), layout='constrained')

axs[0].step(range(151),critical_loads_np_array)
# axs[0].set_xlabel('Tiempo [min]')
axs[0].set_ylabel('Potencia [kW]')
axs[0].set_xticks(np.arange(0, 151, 25))
axs[0].set_yticks(np.arange(0, 4.1, 1))
axs[0].grid(True)

axs[1].step(range(151),HP_loads_np_array)
# axs[1].set_xlabel('Tiempo [min]')
axs[1].set_ylabel('Potencia [kW]')
axs[1].set_xticks(np.arange(0, 151, 25))
axs[1].set_yticks(np.arange(0, 4.1, 1))
axs[1].grid(True)

axs[2].step(range(151),LP_loads_np_array)
axs[2].set_xlabel('Tiempo [min]', fontsize=12)
axs[2].set_ylabel('Potencia [kW]')
axs[2].set_xticks(np.arange(0, 151, 25))
axs[2].set_yticks(np.arange(0, 4.1, 1))
axs[2].grid(True)

plt.show()

def prepare_loads_data(cont, horizon): #EXTRACT LOADS DEMAND FOR SELECTED HORIZON  

    # READ CSV LOADS DEMAND
    critical_loads_dic = pd.read_csv ('D:\\Critical_loads.csv', names = ['critical loads'])
    HP_loads_dic       = pd.read_csv ('D:\\HP_loads.csv', names = ['HP loads'])
    LP_loads_dic       = pd.read_csv ('D:\\LP_loads.csv', names = ['LP loads'])
    
    P_L_crit_dic = {}
    P_L_HP_dic = {}
    P_L_LP_dic = {}
    
    P_L_crit_np_array = critical_loads_dic["critical loads"][cont:(cont+horizon)].to_numpy()
    P_L_HP_np_array   = HP_loads_dic["HP loads"][cont:(cont+horizon)].to_numpy()
    P_L_LP_np_array   = LP_loads_dic["LP loads"][cont:(cont+horizon)].to_numpy()
    
    for k in range(horizon):
        P_L_crit_dic[k] = P_L_crit_np_array[k]
        P_L_HP_dic[k]   = P_L_HP_np_array[k]
        P_L_LP_dic[k]   = P_L_LP_np_array[k]
    
    return P_L_crit_dic, P_L_HP_dic, P_L_LP_dic

#Checks if correct extraction
# (P_L_crit_check, P_L_HP_check, P_L_LP_check) = prepare_loads_data(0, 10)
# print(P_L_crit_check, P_L_HP_check, P_L_LP_check)


########################################## SOLUTIONS ##################################################

def extract_horizon_solution(model):
    
    t      = []
    P_in_aux   = []
    S_L_HP_aux = []
    S_L_LP_aux = []
    SOC_aux    = []
    P_ch_aux   = []
    P_disch_aux= []
    
    for k in model.k:
        t.append(pyo.value(k))
        P_in_aux.append(pyo.value(model.P_in[k]))
        S_L_HP_aux.append(pyo.value(model.S_L[0,k]))
        S_L_LP_aux.append(pyo.value(model.S_L[1,k]))
        SOC_aux.append(pyo.value(model.SOC[k]))
        P_ch_aux.append(pyo.value(model.P_ch[k]))
        P_disch_aux.append(pyo.value(model.P_disch[k]))
    
    return t, P_in_aux, S_L_HP_aux, S_L_LP_aux, SOC_aux, P_ch_aux, P_disch_aux
    

###################################### EXECUTION ########################################

# T = 10
cont = 0
horizon = 150

# P_in    = []
# S_L_HP  = []
# S_L_LP  = []
# SOC     = []
# P_ch    = []
# P_disch = []

# while cont < T:
    
# (P_L_crit_dic, P_L_HP_dic, P_L_LP_dic) = prepare_loads_data(cont, horizon)

# mpc1 = build_model(horizon, P_L_crit_dic, P_L_HP_dic, P_L_LP_dic)

# (t, P_in_aux, S_L_HP_aux, S_L_LP_aux, SOC_aux, P_ch_aux, P_disch_aux) = extract_horizon_solution(mpc1)

#     P_in.append(P_in_aux[0])
#     S_L_HP.append(S_L_HP_aux[0])
#     S_L_LP.append(S_L_LP_aux[0])
#     SOC.append(SOC_aux[0])
#     P_ch.append(P_ch_aux[0])
#     P_disch.append(P_disch_aux[0])

    # cont += 1;
    
# x = np.arange(0, 150)

# P_L_crit = list(P_L_crit_dic.values())
# P_L_HP = list(P_L_HP_dic.values())
# P_L_LP = list(P_L_LP_dic.values())

# plt.figure(figsize=(10, 5)) # POWER BALANCE

# P_generated = []
# for i in range(len(t)): P_generated.append( P_in_aux[i] + P_disch_aux[i] )
# plt.subplot(211) 
# plt.step(x, P_generated, color='green')
# plt.xticks(np.arange(0, 160, 25))
# plt.ylabel('Potencia generada [kW]')
# plt.grid(True)


# P_consumed = []
# for i in range(len(t)): P_consumed.append((-1)*(P_ch_aux[i] + np.multiply(S_L_HP_aux, P_L_HP)[i] + np.multiply(S_L_LP_aux, P_L_LP)[i] + P_L_crit[i]))
# plt.subplot(212) 
# plt.step(x, P_consumed, color='red')
# plt.xlabel('Tiempo [min]')
# plt.xticks(np.arange(0, 160, 25))
# plt.ylabel('Potencia consumida [kW]')
# plt.grid(True)

# plt.show() ##################


# plt.figure(figsize=(15,10)) # POWER DEMAND / SHEDDING

# plt.subplot(311)
# plt.step(x, P_L_crit)

# plt.xticks(np.arange(0, 160, 25))
# plt.xlim(0)
# plt.grid(True)

# plt.subplot(312)
# plt.step(x, P_L_HP)
# plt.step(x, np.multiply(P_L_HP, S_L_HP_aux))
# plt.xticks(np.arange(0, 160, 25))
# plt.xlim(0)
# plt.grid(True)

# plt.subplot(313)
# plt.step(x, P_L_LP)
# plt.step(x, np.multiply(P_L_LP, S_L_LP_aux))
# plt.xticks(np.arange(0, 160, 25))
# plt.xlim(0)
# plt.grid(True)

# plt.show() ######################

# plt.figure(figsize=(10,5)) # SOC
# plt.plot(x, SOC_aux)
# plt.xticks(np.arange(0, 160, 25))
# plt.xlim(0,150)
# plt.ylim(0.2, 1)
# plt.grid(True)

# LO = []
# for i in range(len(x)): LO.append(0.3)
# HI = []
# for i in range(len(x)): HI.append(0.9)
# plt.plot(x, HI, color='red')
# plt.plot(x, LO, color='red')

# plt.show() ######################


# plt.figure(figsize=(10,5)) # CHARGE/DISCHARGE

# plt.step(x, P_ch_aux)
# plt.xticks(np.arange(0, 160, 25))
# plt.xlim(0,150)
# plt.grid(True)

# P_disch_aux_neg = np.multiply(P_disch_aux, -1)
# plt.step(x, P_disch_aux_neg)
# plt.xticks(np.arange(0, 160, 25))
# plt.xlim(0,150)

# plt.show() #######################




############################## PLOT ################################


#plt.figure()

# plt.step(t, S_L_HP_res,'r.-',where='pre')
# # plt.xlabel('Time (hr)')
# # plt.ylabel('Connectors states')
# # plt.xticks(range(0,25,3))
# plt.grid(True)
# plt.show()

    
    
    
    
    
    
    
    
    
    
    
