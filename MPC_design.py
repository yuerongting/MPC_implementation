# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 22:10:14 2021

@author: roy20001
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
import sympy as sym
#from pydmd import DMDc
from scipy.integrate import odeint
import random as rd
import control
from numpy import array, dot
from qpsolvers import solve_qp
# import sympy as sym

# sym.init_printing(use_unicode= True)

interval = 1
t = np.arange(0,43,interval) # days 

def chemo_model(x,t):

    u1 = 1
    u2 = 0.41 
    p1 = 1.25
    p2 = 0.285
    p3 = 1.1
    p4 = 0.12
    p5 = 0.003
    alpha = 0.52
    beta = 0.02
    r = 0.3

    # n = 0
    # tau = 0
    
    # diracDelta = sym.DiracDelta(t - n * tau)
    
    B = x[0]
    E = x[1]
    Ti = x[2]
    Tu = x[3]
    
    
    # u_input = u
    # b=0
    # # u = 0
    # tau = 1 # time interval
    # # if ((t/tau).is_integer() ):
    # if ((t<=1.01)):
    #     b = 1   
    # # else:
    # #     b = 0
    
    # delta_B = b * u_input   #  input

#    dB_dt = -u1 * B - p1 * E * B - p2 * B * Tu
    dB_dt = -u1 * B - p1 * E * B - p2 * B * Tu #+ delta_B
    
    dE_dt = -u2 * E + alpha * Ti + p4 * E * B - p5 * E * Ti
    dTi_dt = -p3 * E * Ti + p2 * B * Tu
    dTu_dt = -p2 * B * Tu + r * (1 - beta * Tu) * Tu
    
#    dE_dt = -u2 * E + alpha * Ti + p4 * E * B - p5 * E * Ti
#    dTi_dt = -p3 * E * Ti + p2 * B * Tu
#    dTu_dt = -p2 * B * Tu + r * (1 - beta * Tu) * Tu
    
    return [dB_dt, dE_dt, dTi_dt, dTu_dt]

def simulate():
    # initial vals
    model_val = np.zeros((len(t),4))
    
    
    B_ini = 0.1
    E_ini= 0.1
    Ti_ini= 0
    Tu_ini= 0.8
    
    # b=1
    b= 0
    
    model_val[0,0] = B_ini
    model_val[0,1] = E_ini
    model_val[0,2] = Ti_ini
    model_val[0,3] = Tu_ini
    model_val_ini = model_val[0]
    
    u=3*np.ones(len(t)*interval)
    
    B_left = []
    drug = []

    for i in range(len(t)-1):
    
        index = int(i * interval)//1
        ts = [t[i],t[i+1]]
        y = odeint(chemo_model, model_val_ini, ts,)
        
        if( (t[i]//1.0).is_integer() ):
            drug_given = b * u[index]
            
            
            y_left_limit = y[1][0]
            B_left.append(y_left_limit)
            
            y[1][0] = y[1][0] + drug_given
            drug.append(drug_given)
            drug.append(0)
            drug.append(0)
        
        
        model_val_ini = y[-1]
        model_val[i+1] = model_val_ini
        # B_left.append(y_left_limit)
    
    return model_val, B_left, drug

def simulate_with_input(initial_val, u):
    # initial vals
    model_val = np.zeros((len(t),4))
    
    
    # b=1
    b= 1
    
    model_val[0,0] = initial_val[0]
    model_val[0,1] = initial_val[1]
    model_val[0,2] = initial_val[2]
    model_val[0,3] = initial_val[3]
    model_val_ini = model_val[0]
    
    u_input = u
    
    B_left = []
    drug = []

    for i in range(len(t)-1):
    
        index = int(i * interval)//1
        ts = [t[i],t[i+1]]
        y = odeint(chemo_model, model_val_ini, ts,)
        
        
        
            
        if( (t[i]//1.0).is_integer() ):
            if(len(u)>i+1):
                
                drug_given = b * u_input[i]
                drug.append(drug_given)
            else:
                drug.append(0)
            
            y_left_limit = y[1][0]
            B_left.append(y_left_limit)
            
            y[1][0] = y[1][0] + drug_given

        
        
        model_val_ini = y[-1]
        model_val[i+1] = model_val_ini
        # B_left.append(y_left_limit)
    
    return model_val, B_left, drug



'''Simulation'''
simulation = []

simulation, B_left, drug = simulate()

initial_val_test = [0.1,0.1,0,0.8]    # test with control input = 1
# control_test = 1
# simulation, B_left, drug = simulate_with_input(initial_val_test, control_test)









'''Test with input'''
# simulation, B_left, drug = simulate_with_input(initial_val_test, u_list)

E_normalize = (simulation[:,1]-min(simulation[:,1]))/(max(simulation[:,1]) - min(simulation[:,1]))
Ti_normalize = (simulation[:,2]-min(simulation[:,2]))/(max(simulation[:,2]) - min(simulation[:,2]))
Tu_normalize = (simulation[:,3]-min(simulation[:,3]))/(max(simulation[:,3]) - min(simulation[:,3]))





'''Original dynamics (normalized cell population)'''
import matplotlib as mat

x_axis = t

# plt.subplot(2,1,1)
# plt.scatter(x_axis, simulation[:,0],label='B^right limit -- vaccine')  

B_left_time = x_axis[1:len(B_left)+1]-0.05
# plt.scatter(B_left_time, B_left, label = 'B^left limit')

B_list = []
B_time = []
for i in range(len(B_left)):
    B_time.append(x_axis[i])
    B_time.append(B_left_time[i])
    
    B_list.append(simulation[:,0][i])
    B_list.append(B_left[i])
    
plt.plot(B_time,B_list, label='B concentration change')
# plt.scatter(x_axis, simulation[:,0],label='B^right limit -- vaccine')  
plt.scatter(B_left_time, B_left, label = 'B^left limit')
plt.legend(loc='upper right', fontsize = 15)




t_interpo = np.zeros(3*(len(t)-1))
for i in range(len(t)-1):
    # if(i >=1):
    t_interpo[3*i] = t[i]
    t_interpo[3*i+1] =  t[i]+0.01
    t_interpo[3*i+2] =  t[i]+0.99
    
drug_act = np.zeros(len(t_interpo) + 3)
for i in range(len(t)-1):    
    drug_act[3*(i+1)] = drug[i]
plt.scatter(  t_interpo  , drug_act[:len(t_interpo)] ,label='u-- Drug administration')
plt.legend(loc='upper right', fontsize = 15)


# plt.subplot(2,1,2)
plt.plot(x_axis, E_normalize, label='E -- Effector cells (10**6)')  
plt.plot(x_axis, Ti_normalize, label='Ti-- Tumor infected cells (10**6)')  
plt.plot(x_axis, Tu_normalize, label='Tu-- Tumor uninfected cells (10**6)')



plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='lower right', fontsize = 15)









''''''
'''Linearization'''
import random as rd
# u1 = 1   # all parameters are perturbed with 10% uncertainty
# u2 = 0.41 
# p1 = 1.25
# p2 = 0.285
# p3 = 1.1
# p4 = 0.12
# p5 = 0.003
# alpha = 0.52
# beta = 0.02
# r = 0.3

import sympy as sym
B, E, Ti, Tu = sym.symbols('B E Ti Tu')
u1, u2, p1, p2, p3, p4, p5, alpha, beta, r = sym.symbols('u1 u2 p1 p2 p3 p4 p5 alpha beta r')


u = 0

dB_dt = -u1 * B - p1 * E * B - p2 * B * Tu + u
dE_dt = -u2 * E + alpha * Ti + p4 * E * B - p5 * E * Ti
dTi_dt = -p3 * E * Ti + p2 * B * Tu
dTu_dt = -p2 * B * Tu + r * (1 - beta * Tu) * Tu

A_jacobian = sym.Matrix([dB_dt,dE_dt,dTi_dt,dTu_dt]).jacobian([B,E,Ti,Tu])

A_jac_num = A_jacobian.subs([('u1', 1), ('u2', 0.41), ('p1', 1.25), ('p2', 0.285), ('p3', 1.1), \
                             ('p4', 0.12), ('p5', 0.003), ('alpha', 0.52), ('beta', 0.02), ('r', 0.3)])





linearized_B = np.array([])
linearized_E = np.array([])
linearized_Ti = np.array([])
linearized_Tu = np.array([])
# A_mat = np.array([])
A_mat = []

B = np.zeros((4,1))
B[0] = 1
# C = np.ones((1,4))
D = 0
sample_int = 5
# period = 50   
for i in range(0,len(t),sample_int):
    day = simulation[i,:]  # simulate model 'period' times
    A = A_jac_num.subs([('B', day[0]), ('E', day[1]), ('Ti', day[2]), ('Tu', day[3])]) # A0, 0; 

    day_initial = day
    # for j in range():
    # j = i    
    t_gap = 1  # simulate every one day
    t_num = 2 # 1 days interval
    
    
    C = [1,0,0,0]
    sys_ss = control.ss(A,B,C,D)
    # t0, res_B = control.initial_response(sys_ss, T = np.linspace(t_gap*(j),t_gap*(j+1),t_num), X0 = day_initial)   
    t0, res_B = control.initial_response(sys_ss, T = np.linspace(0,t_gap,t_num), X0 = day_initial)  
    res_B[res_B<0]=0


    C = [0,1,0,0]
    sys_ss = control.ss(A,B,C,D)
    t0, res_E = control.initial_response(sys_ss, T = np.linspace(0,t_gap,t_num), X0 = day_initial)  
    res_E[0:5] = res_E[0]
    # E_store = res_E[0] * np.ones(5)
    # res_E= E_store
    

    C = [0,0,1,0]
    sys_ss = control.ss(A,B,C,D)
    t0, res_Ti = control.initial_response(sys_ss, T = np.linspace(0,t_gap,t_num), X0 = day_initial)  
    res_Ti[0:5] = res_Ti[0]
    # Ti_store = res_Ti[0] * np.ones(5)
    # res_Ti= Ti_store
    
    
    C = [0,0,0,1]
    sys_ss = control.ss(A,B,C,D)
    t0, res_Tu = control.initial_response(sys_ss, T = np.linspace(0,t_gap,t_num), X0 = day_initial)  
    res_Tu[0:5] = res_Tu[0]
    # Tu_store = res_Tu[0] * np.ones(5)
    # res_Tu= Tu_store
        
    linearized_B = np.concatenate((linearized_B,res_B[0]),axis = None)
    linearized_E = np.concatenate((linearized_E,res_E[0]),axis = None)
    linearized_Ti = np.concatenate((linearized_Ti,res_Ti[0]),axis = None)
    linearized_Tu = np.concatenate((linearized_Tu,res_Tu[0]),axis = None)
    # A_mat = np.concatenate((A_mat, A),axis = None)
    A_mat.append(A)
    
    
    
    
    
    
    







'''Linearization Comparison'''

%matplotlib qt


x_axis = t
x_linear = np.linspace(0,len(x_axis), len(linearized_Tu[0:len(simulation[:,3])])   )
# x_linear = np.linspace(0, len(linearized_B), len(simulation[:,3]))

import matplotlib as mat
test_data = linearized_Tu[0:len(simulation[:,3])]


plt.subplot(2,2,4)
plt.plot(x_axis, simulation[:,3],label='Tu-- original')
# plt.plot(x_axis, dmd.reconstructed_data[3,:].real, '--', label='Tu -- DMD output')
plt.plot(x_linear, linearized_Tu[0:len(simulation[:,3])],label='Tu-- Linearization')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='lower right', fontsize = 15)
# plt.plot(hodmd.dmd_timesteps, hodmd.reconstructed_data[0].real, '--', label='DMD output')

plt.subplot(2,2,3)
plt.plot(x_axis, simulation[:,2],label='Ti-- original')
# plt.plot(x_axis, dmd.reconstructed_data[2,:].real, '--', label='Ti -- DMD output')
plt.plot(x_linear, linearized_Ti[0:len(simulation[:,3])],label='Ti-- Linearization')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='lower right', fontsize = 15)
# plt.plot(hodmd.dmd_timesteps, hodmd.reconstructed_data[0].real, '--', label='DMD output')

plt.subplot(2,2,2)
plt.plot(x_axis, simulation[:,1],label='E-- original')
# plt.plot(x_axis, dmd.reconstructed_data[1,:].real, '--', label='E -- DMD output')
plt.plot(x_linear, linearized_E[0:len(simulation[:,3])],label='E-- Linearization')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='upper right', fontsize = 15)
# plt.plot(hodmd.dmd_timesteps, hodmd.reconstructed_data[0].real, '--', label='DMD output')

plt.subplot(2,2,1)

plt.plot(x_axis, simulation[:,0],label='B-- original')
# plt.plot(x_axis, dmd.reconstructed_data[0,:].real, '--', label='B -- DMD output')
plt.plot(x_linear, linearized_B[0:len(simulation[:,3])],label='B-- Linearization')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='upper right', fontsize = 15)
# plt.plot(hodmd.dmd_timesteps, hodmd.reconstructed_data[0].real, '--', label='DMD output')










'''Parameter setting'''
P = 4
C = np.matrix([1,0,0,1])
D = 0


'''B matrix'''
# B = sym.Matrix([1,0,0,0])    # B matrix
B_mat = []
for i in range(len(x_linear)):
    B_mat.append(B)
    
state_len = len(B)






'''Test'''
# B = sym.Matrix([0,0.0787])
# C = sym.Matrix([[-1,1]])
# R_0 = 0.01
# state_len = 2

# '''Test'''
# # A_fix = A_mat[0]
# # for i in range(len(A_mat)):
# #     A_mat[i] = A_fix
# # B_fix = B_mat[0]
# # for i in range(len(B_mat)):
# #     B_mat[i] = B_fix
# A = np.matrix([[1.1,2],[0,0.95]])
# B = np.matrix([[0],[0.0787]])
# C = np.matrix([-1,1])
# D = 0
# Ts= 1
# A_fix = A
# for i in range(len(A_mat)):
#     A_mat[i] = A_fix
# B_fix = np.matrix([[1],[1]])
# for i in range(len(B_mat)):
#     B_mat[i] = B_fix


'''discrete LQR'''
def dlqr_calculate(G, H, Q, R, returnPE=False):
  '''
  Discrete-time Linear Quadratic Regulator calculation.
  State-feedback control  u[k] = -K*x[k]

  How to apply the function:    
      K = dlqr_calculate(G,H,Q,R)
      K, P, E = dlqr_calculate(G,H,Q,R, return_solution_eigs=True)

  Inputs:
    G, H, Q, R  -> all numpy arrays  (simple float number not allowed)
    returnPE: define as True to return Ricatti solution and final eigenvalues

  Returns:
    K: state feedback gain
    P: Ricatti equation solution
    E: eigenvalues of (G-HK)  (closed loop z-domain poles)
  '''
  from scipy.linalg import solve_discrete_are, inv, eig
  P = solve_discrete_are(G, H, Q, R)  # Ricatti
  K = inv(H.T@P@H + R)@H.T@P@G    #K = (B^T P B + R)^-1 B^T P A 

  if returnPE == False:   return K

  # from numpy.linalg import eigvals
  # eigs = np.array([eigvals(G-H@K)]).T
  return K, P, eigs


'''When t = final time, pick terminal gain K'''
# cost_stage = len(t)-1
cost_stage = 5

# state_space_ini = control.ss(A_mat[cost_stage], B_mat[cost_stage], C, D,Ts)
# (K_lqr, X, E) = control.lqr(state_space_ini, C.T*C, R_0 )

# state_space_ini = control.ss(A,B,C,D,Ts)




'''Q and R matrix'''

Q_0 = C.T * C

# R_0 = np.matrix(0.5)
R_0 = np.matrix(0.5)


cost_stage = len(x_linear)-1
A_type = np.matrix(sym.matrix2numpy(A_mat[cost_stage]))
A = np.matrix(A_type, dtype='float')
# B_type = np.matrix(sym.matrix2numpy(B_mat[cost_stage]))
B_type = np.matrix(B_mat[cost_stage])
B = np.matrix(B_type, dtype='float')


# K_lqr,J,L = control.lqr(A,B,   Q_0,   R_0 )
K_lqr = dlqr_calculate(A,B,   Q_0,   R_0 )

K_ini = np.matrix(K_lqr); 






'''P matrix (Riccati equation, initial value)'''
P_mat = control.dlyap((A_mat[cost_stage] - B * K_ini).T, (Q_0 + K_ini.T * R_0 * K_ini))    # 

# Q = scipy.linalg.block_diag(Q_0,Q_0,Q_0,P_mat)

# R = scipy.linalg.block_diag(R_0,R_0,R_0,R_0)


# Q = scipy.linalg.block_diag(Q_0,Q_0,Q_0,Q_0,Q_0,Q_0,Q_0,Q_0,Q_0,P_mat)

# R = scipy.linalg.block_diag(R_0,R_0,R_0,R_0,R_0,R_0,R_0,R_0,R_0,R_0)

Q = scipy.linalg.block_diag(Q_0,Q_0,Q_0,Q_0,P_mat)

R = scipy.linalg.block_diag(R_0,R_0,R_0,R_0,R_0)





'''Invariant terminal set'''
# import polytope
# X_max = 3
# X_min = 0
# U_max = 5
# U_min = 0
# X = polytope.Polytope(  np.vstack((np.identity(P),np.identity(P)))   , np.array([X_max,X_max,X_max,X_max, -X_min,-X_min,-X_min,-X_min]))
# U = polytope.Polytope(  np.vstack(( np.diag(np.array(K_lqr)[0] ), np.diag(np.array(K_lqr)[0] )    ))   , np.array([U_max,U_max,U_max,U_max   ,   U_min,U_min,U_min,U_min]))


# i=0
# (A_mat[i] + B_mat[i]*K_lqr).inv()

# P = X and U
# while 1 :
#     Pprev = P;
#     P = P and (A_mat[i] + B_mat[i]*K_lqr).inv()*P
#     if Pprev == P:
#         break

# P.plot()





'''Main function'''
u_control_seq = []
simulation_opt = []
B_left_opt = []
drug_opt = []
x_state = simulation
ff = 0
P = 5

initial_val = np.array([0.1,0.1,0,0.8])
x_0 = np.matrix(initial_val).T
u_list = []
sys_res = []
sys_res.append(x_0)



# for ff in range(len(t)-P-1):

# for ff in range(20):
    
    ff= 0
    
    
    '''M matrix'''
    M_list = []
    # M = sym.Matrix([])
    M = np.matrix([])
    
    
    for i in range(0+ff,P+ff):
        A_matrixA = 1
        for j in range(i+1):
            # print(j)
            A_matrixA =  A_mat[j] * A_matrixA   # A[0] ~ A[19]   descending:  A[19] * A[18] * ...
        M_list.append(A_matrixA) 
    
    
    # M = sym.Matrix([M_list[0],M_list[1]])
    # M = np.array(M_list[0]).astype(np.float64)
    M = np.array(M_list).reshape(-1,state_len)

    '''Cc matrix'''
    
    # P =20
    Cc = sym.zeros(rows = state_len*P,cols = P)
    # Cc = sym.Matrix([])
    
    AA_list=[1]
    # ff = 0
    for n in range(1+ff, P + ff):

        AA_list.append(A_mat[n])     # A[1]~ [19]   [2~20]  [3~ 21]
        # AA_list.append(n) 
    
    
    # P = 3
    
    for i in range(P):   # row        
        
        current_line = [] 
        current_element = 1
        for j in range(i):  # col
            # start_element = AA_list[j+1]
            
            # print(i-j)
            current_element =  current_element * AA_list[i-j]
            # current_element = start_element * AA_list[i-j]
            current_line.append(current_element)
            
        current_line = current_line[len(current_line)::-1]
            

                
        for col in range(i):
            # print(col)

            Cc[ state_len*i : state_len*(i+1) , col ] = current_line[col] * B_mat[0 + ff]    # times B matrix
        # for j in range(i+1):
            
        Cc[ state_len*i : state_len*(i+1) , i ] = B_mat[0 + ff]
        
 
    '''H and F and G matrix'''
    H = Cc.T * Q * Cc + R   # calculate H , F, G for each sampling time
    F = Cc.T * Q * M
    G = M.T @ Q @ M + Q_0
    
    
    '''Global optima control input (H matrix)'''
    def is_pos_def(x):
        return np.all(np.linalg.eigvals(x) > 0)
    
    
    H.eigenvals()
    
    
    '''Chose 1st resulting gain'''
    gain_U = H.inv()*F 
    gain_current = gain_U[0,:]
    
    '''
    Test
    '''
    # initial_val_test = np.array([0.5,-0.5])
    
    
    
    initial_val = np.array([0.1,0.1,0,0.8])
    # u_optimal = -(gain_current * (np.matrix(initial_val).T))[0]   # u[1] = K * x[0]
    u_optimal = -(gain_current * (np.matrix(initial_val_test).T))[0]   # u[1] = K * x[0]
    u_control_seq.append(u_optimal)   # current control input
    
    
    
    
    
    '''
    system test
    '''
    # Ts = 1
    # C_observe_all = np.eye(state_len)
    # new_ss = control.ss(A_mat[ff],B_mat[ff],C_observe_all,0,Ts)
    # new_sys = C*control.feedback(new_ss,gain_current)
    
    
    # x_0 = np.matrix(initial_val_test).T  # update current initial state
    
    # end = P
    # time_simu = np.linspace(0, end, end+1)
    # u_unconstrained = np.zeros(len(time_simu))
    # import control.matlab
    
    
    # y_out, Time, xout = control.matlab.lsim(new_sys, u_unconstrained,time_simu, X0 = initial_val_test)
    # # x_current_state = y_out
    # # plt.scatter(Time, Output)
    # # plt.scatter(Time, xout[:,0], label = 'B')
    # # plt.plot(Time, xout[:,1], label = 'E')
    # # plt.plot(Time, xout[:,2], label = 'Ti')
    # # plt.plot(Time, xout[:,3], label = 'Tu')
    # # plt.plot(Time, y_out, label = 'Output')
    # # plt.legend()
    
    # '''Test output'''
    # plt.plot(Time, xout[:,0], label = 'x1')
    # plt.plot(Time, xout[:,1], label = 'x2')
    
    # plt.plot(Time, y_out, label = 'Output')
    # plt.legend()
    
    
    # xout[:,0][1]  # index [1] current state
    # xout[:,1][1]
    
    
    
    '''
    '''
    
    
    # initial_val = x_state[0,:]
    
    # simulation_current, B_left_current, drug_current = simulate_with_input(initial_val, u_optimal)   # not ode45, use state space
    
    # simulation_opt.append(simulation_current[1,:])   # output x[1]
    # B_left_opt.append(B_left_current[0])
    # drug_opt.append(drug_current)
    
    
    # x_state = np.zeros([P, len(simulation[0,:])] )
    # x_state[0,:] = simulation_current[1,:]


    
    Ac_state = np.zeros([2*(Cc.shape[0]), Cc.shape[1]])
    bc_state = np.zeros([(2*Cc.shape[0]), 1])
    # Bc_xk = np.zeros([(2*Cc.shape[0]), A_mat[0].shape[1]])
    Bc_xk = np.zeros([2*Cc.shape[0],1])
    
    
    
    
    
    import qpsolvers
    
    " Cholesky decomposition "
    P_para_mat = np.array(H, dtype='float')
    # P_chol = np.linalg.cholesky(P_para_mat)
    # P_para = np.linalg.inv(np.matrix(P_chol))*np.eye(H.shape[0])

    q_para = np.array(2* F * x_0, dtype='float').T[0]
    
    
    
    
    
    
    '''Constraint setting'''
    state_max = 4*np.matrix(np.ones(P)).T
    state_min = -0.3*np.matrix(np.ones(P)).T
    input_max = 5.9*np.matrix(np.ones(P)).T
    input_min = -1*np.matrix(np.ones(P)).T
    # state_constraint_index = 2*state_len
    
    x_0 = np.matrix(initial_val_test).T 
    Ac_state = np.zeros([2*(P+len(B)), P])  #  P + P + len(B) + len(B) 
    bc_state = np.zeros([2*(P+len(B)), 1])
    # Bc_xk = np.zeros([(2*Cc.shape[0]), A_mat[0].shape[1]])
    Bc_xk = np.zeros([2*(P+len(B)),1])    
    
    
    
    
    # x_0 = np.matrix(initial_val).T
    u_list = []
    sys_res = []
    # sys_res.append(x_0)
    
    
    for i in range(P):
        
        
        # i=0
        
        
        Ac_state[ 0:P , :] = np.eye(P)
        Ac_state[ P:2*P, :] = - np.eye(P)
        
        Ac_state[ 2*P : 2*P + state_len , :] = Cc[ state_len*i : state_len*(i+1), :]
        Ac_state[ 2*P + state_len : 2*P + 2*state_len, :] =-Cc[ state_len*i : state_len*(i+1), :]
        
        
        # bc_state[ 4*i : 4*i+state_len] = state_max[i]
        # bc_state[ 4*i+state_len : 4*(i+1)] = -state_max[i]
        bc_state[ 0:P , :] = input_max[i]
        bc_state[ P:2*P, :] = -input_min[i]
        
        bc_state[ 2*P : state_len , :] = state_max[i]
        bc_state[ 2*P+state_len : 2*P+2*state_len, :] = -state_min[i]
        
        
        Bc_xk[ 0:P ,:] = np.zeros([P,1])
        Bc_xk[ P:2*P, :] = np.zeros([P,1])
        
        Bc_xk[ 2*P : 2*P+state_len , :] = - M[ state_len * i : state_len*(i+1),:] * x_0
        Bc_xk[ 2*P+state_len : 2*P+2*state_len, :] = + M[ state_len * i : state_len*(i+1),:] * x_0
        
        # Bc_xk[ 4*i : 4*i+state_len] = M[2*i:2*(i+1),:] * x_0
        # Bc_xk[ 4*i+state_len : 4*(i+1)] = -M[2*i:2*(i+1),:] * x_0
        b_constraint = (bc_state+Bc_xk).reshape(2*P + 2*state_len,)
        
        
        # index_from = i*2
        # index_end = (i+1)*2
        # u_opt_constrained = qpsolvers.solve_qp(P_para_mat,q_para,Ac_state[index_from:index_end,:],b_constraint[index_from:index_end])
        # u_opt_constrained = qpsolvers.solve_qp(P_para_mat,q_para, Ac_state[0:8],b_constraint[0:8],initvals = x_0)
        
        
        u_opt_constrained = qpsolvers.solve_qp(P_para_mat,q_para, Ac_state,b_constraint,initvals = x_0)
        # u_opt_constrained
        
        
        # u_opt_constrained = qpsolvers.osqp_solve_qp(P_para_mat,q_para, Ac_state,b_constraint,initvals = x_0, eps_abs=0.01, eps_rel=0.01)
        # u_opt_constrained = qpsolvers.osqp_solve_qp(P_para_mat,q_para)
        
        import cvxopt.solvers
        from cvxopt import matrix
        cvxopt.solvers.options['maxiters'] = 200
        cvxopt.solvers.options['feastol'] = 1e-3
        
        result = cvxopt.solvers.qp(matrix(P_para_mat), matrix(q_para), matrix(Ac_state), matrix(b_constraint))
        result['x'][0]
        result['x'][1]
        result['x'][2]
        result['x'][3]
        result['x'][4]
        
        # qpsolvers.scs_solve_qp(P_para_mat,q_para)
        # u_opt_constrained = qpsolvers.solve_qp(P_para_mat,q_para,Ac_state,b_constraint)
        
        
        C_observe_all = np.eye(state_len)
        sys_test = C_observe_all * (A_mat[i+ff] * x_0 + B_mat[i+ff] * result['x'][0])
        # sys_test = C_observe_all * (A_mat[i+ff] * x_0 + B_mat[i+ff] * u_opt_constrained[0])
        # sys_test = (  A_mat[i+ff] * x_0 + B_mat[i+ff] * u_opt_constrained[0])
        sys_test_out = C_observe_all * sys_test
        # u_opt_constrained[0]
        
        
        x_0 = sys_test

        u_list.append(result['x'][0])
        # u_list.append(u_opt_constrained[0])
        sys_res.append(sys_test)        
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
   
    
    # from cvxopt import matrix, solvers
    

    
    # " Cholesky decomposition "
    # P_para_mat = np.array(H, dtype='float')
    # P_chol = np.linalg.cholesky(P_para_mat)
    # P_para = np.linalg.inv(np.matrix(P_chol))*np.eye(H.shape[0])



    # q_para = np.array(2* F * x_0, dtype='float').T[0]
    
    
    # # Ac = np.zeros([Cc.shape[], len(Cc)])
    
    # import qpsolvers
    # import quadprog
    # u_opt_constrained = quadprog(P_para,q_para)
    # u_opt_constrained = qpsolvers.solve_qp(P_para,q_para,)
    
    
    # # G_para = Ac_state
    # # h_para = b0_state    
    # Ac_list = []
    # for i in range(1,int(Cc.shape[0]/P)+1):
    #     # b0_state[2*i:2*i+2] = state_max - M[i]*x_0
    #     # b0_state[i*2+1] = -state_max + M[i]*x_0
    #     state_len = 2
        
    #     b0_max = state_max - M[2*(i-1):2*i]*x_0
    #     b0_min = -state_max + M[2*(i-1):2*i]*x_0
    #     b0 = np.array(np.vstack((b0_max,b0_min))[:,0]).reshape(-1)
        
        



        
    #     Ac_pos = np.matrix(Cc[2*(i-1):2*i,:],dtype='float')
    #     Ac_neg = np.matrix(- Cc[2*(i-1):2*i,:],dtype='float')
    #     Ac = np.array(np.vstack((Ac_pos,Ac_neg)))
    #     # Ac_neg = np.array(Ac_neg)
    #     Ac_list.append(Ac)
        
        
    #     ''' test '''
    #     Ac = np.zeros([8,4])       
    #     b0 = np.zeros([8,1])

        
        
        
        
        
    #     # u_opt_constrained = qpsolvers.solve_qp(P_para,q_para,Ac,b0)
    #     # u_opt_constrained = qpsolvers.solve_qp(P_para,q_para, G = Ac, h = b0, lb=lb_para, ub=ub_para)
        
        
    #     # print(qpsolvers.available_solvers)

    # # A_para = 0
    # # b_para = 0
    # lb_para = -10*np.ones(int(P*2))
    # ub_para = 10*np.ones(int(P*2))
    # # u_opt_constrained = qpsolvers.solve_qp(P_para,q_para,G_para,h_para)
    
    
    
    # qpsolvers.osqp_solve_qp(Ac,b0)
    # # u_opt_constrained = qpsolvers.solve_qp(0.5*P_para, 0.5*q_para, Ac, b0, lb=lb_para, ub=ub_para  )####################
    
    
    # def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    #     qp_G = .5 * (P + P.T)   # make sure P is symmetric
    #     qp_a = -q
    #     if A is not None:
    #         qp_C = -numpy.vstack([A, G]).T
    #         qp_b = -numpy.hstack([b, h])
    #         meq = A.shape[0]
    #     else:  # no equality constraint
    #         qp_C = -G.T
    #         qp_b = -h
    #         meq = 0
    #     return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]
    
    # u_opt_constrained = quadprog_solve_qp(P_para_mat.astype('double'),q_para.astype('double'),Ac.astype('double'),b0.astype('double'))
    # # u_opt_constrained = quadprog.solve_qp(P_para_mat.astype('double'),q_para.astype('double'),Ac.astype('double'),b0.astype('double'))

    




    # M = array([[1., 2., 0.], [-8., 3., 2.], [0., 1., 1.]])
    # P = dot(M.T, M)  # this is a positive definite matrix
    # q = dot(array([3., 2., 3.]), M).reshape((3,))
    # G = array([[1., 2., 1.], [2., 0., 1.], [-1., 2., -1.]])
    # h = array([3., 2., -2.]).reshape((3,))
    # # A = array([1., 1., 1.])
    # # b = array([1.])
    # u_opt_constrained = qpsolvers.solve_qp(P,q,G,h)
    
    

    
'''Optimized dynamics (normalized cell population)'''
B_right_opt = np.array(simulation_opt)[:,0]
x_axis_opt = x_axis[:len(B_right_opt)]
plt.subplot(2,1,1)
plt.scatter(x_axis_opt, B_right_opt, label='B^right limit -- vaccine')  
plt.scatter(x_axis_opt, B_left_opt, label = 'B^left limit')



t_interpo_opt = np.zeros(3*(len(x_axis_opt)-1))
for i in range(len(x_axis_opt)-1):
    # if(i >=1):
    # t_interpo_opt[3*i] = t[i+len(t)-P-1]
    t_interpo_opt[3*i] = t[i]
    t_interpo_opt[3*i+1] =  t[i]+0.01
    t_interpo_opt[3*i+2] =  t[i]+0.99
    
drug_act_opt = np.zeros(len(t_interpo_opt) + 3)

for i in range(len(x_axis_opt)-1):
    drug_act_opt[3+3*i] = u_control_seq[i]
plt.plot(  t_interpo_opt  , drug_act_opt[:len(t_interpo_opt)] ,label='u-- Drug administration')
plt.legend(loc='upper right', fontsize = 15)







plt.subplot(2,1,2)
simulation_opt_array = np.array(simulation_opt)

E_nor_opt = simulation_opt_array[:,1] / (max(simulation_opt_array[:,1])- min(simulation_opt_array[:,1]))
Ti_nor_opt = simulation_opt_array[:,2] / (max(simulation_opt_array[:,2])- min(simulation_opt_array[:,2]))
Tu_nor_opt = simulation_opt_array[:,3] / (max(simulation_opt_array[:,3])- min(simulation_opt_array[:,3]))

plt.plot(x_axis_opt, E_nor_opt, label='E -- Effector cells (10**6)')  
plt.plot(x_axis_opt, Ti_nor_opt, label='Ti-- Tumor infected cells (10**6)')  
plt.plot(x_axis_opt, Tu_nor_opt, label='Tu-- Tumor uninfected cells (10**6)')



plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='upper right', fontsize = 15)











'''Reference trajectory of Tu'''
B_ini = 0.1
E_ini= 0.1
Ti_ini= 0
Tu_ini= 0.8
tau = 10

Tu_traj = Tu_ini * np.exp(-t/tau)

plt.plot(t, Tu_traj, label = 'Tu Reference (exp decay)')
plt.legend()




'''Cost func'''



# i = 0
# J_k = x_state[0].T * Q_0 * x_state[0] + x_state[0:P].T * Q * x_state[0:P] + U.T * R * U




'''Simulation'''

x_state = np.zeros([P, len(simulation[0,:])] )
x_state[:,0] = simulation[:P,0]
x_state[:,1] = simulation[:P,1]
x_state[:,2] = simulation[:P,2]
x_state[:,3] = simulation[:P,3]

simulation = []
PQ_update  = simulate()
simulation = PQ_update
simulation[:,0] = simulation[:,0]/(max(simulation[:,0]) - min(simulation[:,0]))
simulation[:,1] = simulation[:,1]/(max(simulation[:,1]) - min(simulation[:,1]))
simulation[:,2] = simulation[:,2]/(max(simulation[:,2]) - min(simulation[:,2]))
simulation[:,3] = simulation[:,3]/(max(simulation[:,3]) - min(simulation[:,3]))
simulation[:,4] = simulation[:,4]


'''Original dynamics (normalized)'''
import matplotlib as mat

x_axis = t

plt.plot(x_axis, simulation[:,0],label='B -- vaccine (doses)')  
plt.plot(x_axis, simulation[:,1],label='E -- Effector cells (10**6)')  
plt.plot(x_axis, simulation[:,2],label='Ti-- Tumor infected cells (10**6)')  
plt.plot(x_axis, simulation[:,3],label='Tu-- Tumor uninfected cells (10**6)'),  
plt.plot(x_axis, simulation[:,4],label='Drug dose (3~7)'),  
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='upper right', fontsize = 15)










from pydmd import DMD,DMDc
import scipy.integrate
'''DMD'''
X = simulation
fig = plt.plot(scipy.linalg.svdvals(np.array([X.flatten() for X in X]).T), 'o')
# fig = plt.plot(scipy.linalg.svdvals(np.array([simulation.flatten() for simulation in simulation]).T), '*')

from pydmd import MrDMD,OptDMD
# dmd = MrDMD(DMD())

dmd = DMD(svd_rank=3, tlsq_rank=0, exact=True, opt=True)  # not truncated, rank truncation of total least square = 2
# dmd = OptDMD(factorization='evd', svd_rank=0, tlsq_rank=0, opt=False)


dmd.fit(simulation)
dmd.plot_eigs()

# dmd_states = [state.reshape(simulation[:,0].shape) for state in dmd.reconstructed_data.T]

Atilde = dmd.atilde

mode = dmd.modes


# reconstructed_data_0 = abs(dmd.reconstructed_data[:,0])
# reconstructed_data_1 = abs(dmd.reconstructed_data[:,1])
# reconstructed_data_1[reconstructed_data_1<0] = 0
# reconstructed_data_2 = abs(dmd.reconstructed_data[:,2])
# reconstructed_data_2[reconstructed_data_2<0] = 0
# reconstructed_data_3 = abs(dmd.reconstructed_data[:,3])
# reconstructed_data_3[reconstructed_data_3<0] = 0
# reconstructed_data_0 = abs(dmd.reconstructed_data[:,0])


reconstructed_data_0 = dmd.reconstructed_data[:,0].real
# reconstructed_data_0[reconstructed_data_0<0] = 0
reconstructed_data_1 = dmd.reconstructed_data[:,1].real
reconstructed_data_1[reconstructed_data_1<0] = 0
reconstructed_data_2 = dmd.reconstructed_data[:,2].real
reconstructed_data_2[reconstructed_data_2<0] = 0
reconstructed_data_3 = dmd.reconstructed_data[:,3].real
reconstructed_data_3[reconstructed_data_3<0] = 0

'''DMD comparison'''
x_axis = t
plt.subplot(2,2,4)
plt.plot(x_axis, simulation[:,3],label='Tu-- original')
plt.plot(x_axis, reconstructed_data_3, '--', label='Tu -- DMD output')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='lower right', fontsize = 15)

plt.subplot(2,2,3)
plt.plot(x_axis, simulation[:,2],label='Ti-- original')
plt.plot(x_axis, reconstructed_data_2, '--', label='Ti -- DMD output')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='lower right', fontsize = 15)

plt.subplot(2,2,2)
plt.plot(x_axis, simulation[:,1],label='E-- original')
plt.plot(x_axis, reconstructed_data_1, '--', label='E -- DMD output')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='upper right', fontsize = 15)

plt.subplot(2,2,1)
plt.plot(x_axis, simulation[:,0],label='B-- original')
plt.plot(x_axis, reconstructed_data_0, '--', label='B -- DMD output')
plt.xlabel('Time (day)', fontproperties=mat.font_manager.FontProperties(size=15))
plt.ylabel('Normailized range', fontproperties=mat.font_manager.FontProperties(size=15))
plt.legend(loc='upper right', fontsize = 15)












'''Drug schedule'''


    
    

























