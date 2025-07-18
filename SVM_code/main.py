X = np.array([
    [27.65, 15.65],
    [23.1 , 14.6 ],
    [23.5 , 15.2 ],
    [24.05, 14.9 ],
    [24.5 , 14.7 ],
    [14.15, 17.35],
    [14.3 , 16.8 ],
    [14.3 , 15.75],
    [14.75, 15.1 ],
    [15.35, 15.5 ]
])

y = np.array([1, 1, 1, 1, 1, -1, -1, -1, -1, -1])


epochs = 500
lr = 0.0001 # lr * epochs >=1
time_step = 100
ct = 1/500


# %%
def dual_rail_encoding_arr(x):
    x_p = []
    x_n = []
    for i in x:
        p = []
        n = []

        for j in i:
            if j > 0:
                p.append(j)
                n.append(0)
            else:
                p.append(0)
                n.append(-j)
        x_p.append(p)
        x_n.append(n)
    return x_p, x_n

def dual_rail_encoding_num(x):
    x_p = []
    x_n = []

    for i in x:
        if i > 0:
            x_p.append(i)
            x_n.append(0)
        else:
            x_p.append(0)
            x_n.append(-i)
    return x_p, x_n


# %%
X_p, X_n = dual_rail_encoding_arr(X)
Y_p, Y_n = dual_rail_encoding_num(y)
print(X_p, X_n)
print(Y_p, Y_n)

# %%
from scipy.integrate import odeint

# %%
def osciallation(num_species, time_Duration = 165):
    O_t = np.zeros(num_species)
    O_t[0] = 1
    t = np.linspace(0, time_Duration, time_Duration+1)
    dope = 10 ** (-4)
    d_o_T = np.ones(num_species) * dope

    def doped_model(O_t, t,  d_o_T):
        dO_t = np.zeros(num_species)
        dO_t[0] = -O_t[0]*O_t[1] + O_t[num_species-1]*O_t[0] + d_o_T[0]* O_t[num_species-1] - d_o_T[1]*O_t[0]
        for i in range(1, num_species-1):
            dO_t[i] = O_t[i-1]*O_t[i] - O_t[i]*O_t[i+1] + d_o_T[i]*O_t[i-1] - d_o_T[i+1]*O_t[i]
        dO_t[num_species-1] = O_t[num_species-2]*O_t[num_species-1] - O_t[num_species-1]*O_t[0] + d_o_T[num_species-1]*O_t[num_species-2] - d_o_T[0]*O_t[num_species-1]
        return dO_t
    
    O = odeint(doped_model, O_t, t, args=(d_o_T,))
    return O


sol = osciallation(18)

# plotting the solution 
import matplotlib.pyplot as plt
plt.plot(sol)
# plt.legend([f"O{i}" for i in range(1, 19)])
plt.show()

# %%
# laoding function
def load_f(y,t,a):
    c = y[0]
    dc_dt = -c + a
    return [ dc_dt]

def load_wrapper(a, c, time_start):
    y0 = [c]
    # print(a)
    t = np.linspace(time_start, time_start+time_step, time_step+1)
    
    sol = odeint(load_f, y0, t, args=(a,))

    # a_f = sol[:, 0]
    # b_f = sol[:, 1]
    c_f = sol[:, 0]
    return c_f[-1]


# %%
# multiplication function
def product_f(y,t,a,b):
    c = y[0]
    dc_dt = -c + a*b
    return [ dc_dt]

def product_wrapper(a,b,c, time_start):
    y0 = [c]
    t = np.linspace(time_start, time_start+time_step, time_step+1)
    
    sol = odeint(product_f, y0, t, args=(a,b))

    # a_f = sol[:, 0]
    # b_f = sol[:, 1]
    c_f = sol[:, 0]
    return c_f[-1]

# %%
# sum
def sum_f(y,t,a,b):
    c = y[0]
    dc_dt = -c + a + b
    return [ dc_dt]

def sum_wrapper(a,b,c, time_start):
    y0 = [c]
    t = np.linspace(time_start, time_start+time_step, time_step+1)
    
    sol = odeint(sum_f, y0, t, args=(a,b))

    # a_f = sol[:, 0]
    # b_f = sol[:, 1]
    c_f = sol[:, 0]
    return c_f[-1]

# %%
# division
def div_f(y,t,a,b):
    c = y[0]
    dc_dt = a - b*c
    return [ dc_dt]

def div_wrapper(a,b,c, time_start):
    y0 = [c]
    t = np.linspace(time_start, time_start+time_step, time_step+1)
    
    sol = odeint(div_f, y0, t, args=(a,b))

    # a_f = sol[:, 0]
    # b_f = sol[:, 1]
    c_f = sol[:, 0]
    return c_f[-1]


# %%
# comparison
def comparison_f(y, t, a, b):
    dydt = []
    bgta = y[0]
    blta = y[1]

    temp = b*blta - bgta*a
    dydt.append(temp)

    temp = bgta*a - b*blta
    dydt.append(temp)

    return dydt

def Approx_majority (y, t):
    kgtq = y[0]
    kltq = y[1]
    b_help = y[2]

    dydt = []
    
    temp =  - kgtq*kltq + kgtq *b_help
    dydt.append(temp)

    temp = kltq*b_help - kgtq*kltq
    dydt.append(temp)

    temp = 2* kgtq*kltq - kgtq*b_help - kltq*b_help 
    dydt.append(temp)

    return dydt



def comparison_wrapper(a,b,c, time_start):

    bgta = c
    blta = c

    y0 = [bgta, blta]
    t = np.linspace(time_start, time_start+time_step, time_step+1)
    
    sol = odeint(comparison_f, y0, t, args=(a,b))

    # a_f = sol[:, 0]
    # b_f = sol[:, 1]
    bgta = sol[:, 0][-1]
    blta = sol[:, 1][-1]
    help_b = 0

    # print(bgta, blta, help_b)



    y_0_1 = [bgta, blta, help_b]

    t = np.linspace(time_start + time_step , time_start + 2*time_step, time_step+1)

    sol = odeint(Approx_majority, y_0_1, t)
    bgta = sol[:,0][-1]
    blta = sol[:,1][-1]
    help_b = sol[:,2][-1]

    return bgta


# %%
import numpy as np
from sklearn import svm




# %%
# weight species W_1 , W_2 (set by looking at the num of features)
# bias species B

w_0_p = 1
w_0_n = 0
w_1_p = 1
w_1_n = 0
b_p = 0
b_n = 0
P_p = 0
P_n = 0
Q_p = 0
Q_n =0

# intermediates
w_0_p_x_p_0 = 0
w_1_p_x_p_1 = 0
w_0_n_x_n_0 = 0
w_1_n_x_n_1 = 0
w_0_p_x_n_0 = 0
w_0_n_x_p_0 = 0
w_1_p_x_n_1 = 0
w_1_n_x_p_1 = 0

# storing the values of intermediates
w_0_p_x_p_0_arr = [w_0_p_x_p_0]
w_1_p_x_p_1_arr = [w_1_p_x_p_1]
w_0_n_x_n_0_arr = [w_0_n_x_n_0]
w_1_n_x_n_1_arr = [w_1_n_x_n_1]
w_0_p_x_n_0_arr = [w_0_p_x_n_0]
w_0_n_x_p_0_arr = [w_0_n_x_p_0]
w_1_p_x_n_1_arr = [w_1_p_x_n_1]
w_1_n_x_p_1_arr = [w_1_n_x_p_1]



w_0_p_arr = [w_0_p]
w_0_n_arr = [w_0_n]
w_1_p_arr = [w_1_p]
w_1_n_arr = [w_1_n]
b_p_arr = [b_p]
b_n_arr = [b_n]
P_p_arr = [P_p]
P_n_arr = [P_n]
Q_p_arr = [Q_p]
Q_n_arr = [Q_n]
loss_arr = [0]



time_index =0 

for t in range(epochs):
    time_start = t*time_step

    for i, x in enumerate(X):
        # defining species here as complexes
        x_p = X_p[i]
        x_p_0 = x_p[0]
        x_p_1 = x_p[1]

        x_n = X_n[i]
        x_n_0 = x_n[0]
        x_n_1 = x_n[1]

        y_p = Y_p[i]
        y_n = Y_n[i]

        x_p_0_y_p = product_wrapper(x_p_0, y_p, 0, time_start)
        x_p_1_y_p = product_wrapper(x_p_1, y_p, 0, time_start)
        x_n_0_y_n = product_wrapper(x_n_0, y_n, 0, time_start)
        x_n_1_y_n = product_wrapper(x_n_1, y_n, 0, time_start)
        x_p_0_y_n = product_wrapper(x_p_0, y_n, 0, time_start)
        x_p_1_y_n = product_wrapper(x_p_1, y_n, 0, time_start)
        x_n_0_y_p = product_wrapper(x_n_0, y_p, 0, time_start)
        x_n_1_y_p = product_wrapper(x_n_1, y_p, 0, time_start)

        # P_p = w_0_p * x_p[0] + w_1_p * x_p[1] + w_0_n * x_n[0] + w_1_n * x_n[1] + b_n
        w_0_p_x_p_0 = product_wrapper(w_0_p, x_p_0, w_0_p_x_p_0, time_start)
        w_1_p_x_p_1 = product_wrapper(w_1_p, x_p_1, w_1_p_x_p_1, time_start)
        w_0_n_x_n_0 = product_wrapper(w_0_n, x_n_0, w_0_n_x_n_0, time_start)
        w_1_n_x_n_1 = product_wrapper(w_1_n, x_n_1, w_1_n_x_n_1, time_start)
        P_p = sum_wrapper(w_0_p_x_p_0, w_1_p_x_p_1, P_p, time_start)
        P_p = sum_wrapper(P_p, w_0_n_x_n_0, P_p, time_start)
        P_p = sum_wrapper(P_p, w_1_n_x_n_1, P_p, time_start)
        P_p = sum_wrapper(P_p, b_n, P_p, time_start)
        
        # P_n = w_0_p * x_n[0] + w_1_p * x_n[1] + w_0_n * x_p[0] + w_1_n * x_p[1] + b_p
        w_0_p_x_n_0 = product_wrapper(w_0_p, x_n_0, w_0_p_x_n_0, time_start)
        w_0_n_x_p_0 = product_wrapper(w_0_n, x_p_0, w_0_n_x_p_0, time_start)
        w_1_p_x_n_1 = product_wrapper(w_1_p, x_n_1, w_1_p_x_n_1, time_start)
        w_1_n_x_p_1 = product_wrapper(w_1_n, x_p_1, w_1_n_x_p_1, time_start)
        P_n = sum_wrapper(w_0_p_x_n_0, w_1_p_x_n_1, P_n, time_start)
        P_n = sum_wrapper(P_n, w_0_n_x_p_0, P_n, time_start)
        P_n = sum_wrapper(P_n, w_1_n_x_p_1, P_n, time_start)
        P_n = sum_wrapper(P_n, b_p, P_n, time_start)
        
        temp_1 = 1
        temp_1 = product_wrapper(y_p, P_p, temp_1, time_start)

        temp_2 = 1
        temp_2 = product_wrapper(y_n, P_n, temp_2, time_start)

        temp_3 = 1
        temp_3 = product_wrapper(y_p, P_n, temp_3, time_start)

        temp_4 = 1
        temp_4 = product_wrapper(y_n, P_p, temp_4, time_start)

        # Q_p = Y_p[i]* P_p + Y_n[i] * P_n
        Q_p = sum_wrapper(temp_1, temp_2, Q_p, time_start)
        # Q_n = Y_p[i]* P_n + Y_n[i] * P_p
        Q_n = sum_wrapper(temp_3, temp_4, Q_n, time_start)

        Q_n_t = 1
        Q_n_t = sum_wrapper(Q_n, 1, Q_n_t, time_start)

        const = lr *ct
        const_t = 0

        bgta = comparison_wrapper(Q_p, Q_n_t,0.5,time_start)

        # w_0_p = w_0_p + lr * (2 * ct * w_0_n) + lr * (x_p[0] * Y_p[i] + x_n[0] * Y_n[i]) * bgta
        const_1 = product_wrapper(w_0_n, const, const_t, time_start)
        w_0_p = sum_wrapper(w_0_p, const_1, w_0_p, time_start)
        const_2 = sum_wrapper(x_p_0_y_p, x_n_0_y_n, const_t, time_start)
        const_2 =product_wrapper(const_2, lr, const_t, time_start)
        const_2 = product_wrapper(const_2, bgta, const_2, time_start)
        w_0_p = sum_wrapper(w_0_p, const_2, w_0_p, time_start)

        # w_1_p = w_1_p + lr * (2 * ct * w_1_n) + lr * (x_p[1] * Y_p[i] + x_n[1] * Y_n[i]) * bgta
        const_3 = product_wrapper(w_1_n, const, const_t, time_start)
        w_1_p = sum_wrapper(w_1_p, const_3, w_1_p, time_start)  
        const_4 = sum_wrapper(x_p_1_y_p, x_n_1_y_n, const_t, time_start)
        const_4 =product_wrapper(const_4, lr, const_t, time_start)
        const_4 = product_wrapper(const_4, bgta, const_4, time_start)
        w_1_p = sum_wrapper(w_1_p, const_4, w_1_p, time_start)

        # w_0_n = w_0_n + lr * (2 * ct * w_0_p) + lr * (x_n[0] * Y_p[i] + x_p[0] * Y_n[i]) * bgta
        const_5 = product_wrapper(w_0_p, const, const_t, time_start)
        w_0_n = sum_wrapper(w_0_n, const_5, w_0_n, time_start)
        const_6 = sum_wrapper(x_n_0_y_p, x_p_0_y_n, const_t, time_start)
        const_6 =product_wrapper(const_6, lr, const_t, time_start)
        const_6 = product_wrapper(const_6, bgta, const_6, time_start)
        w_0_n = sum_wrapper(w_0_n, const_6, w_0_n, time_start)

        # w_1_n = w_1_n + lr * (2 * ct * w_1_p) + lr * (x_n[1] * Y_p[i] + x_p[1] * Y_n[i]) * bgta
        const_7 = product_wrapper(w_1_p, const, const_t, time_start)
        w_1_n = sum_wrapper(w_1_n, const_7, w_1_n, time_start)
        const_8 = sum_wrapper(x_n_1_y_p, x_p_1_y_n, const_t, time_start)
        const_8 =product_wrapper(const_8, lr, const_t, time_start)
        const_8 = product_wrapper(const_8, bgta, const_8, time_start)
        w_1_n = sum_wrapper(w_1_n, const_8, w_1_n, time_start)

        # b_p = b_p + lr * Y_n[i] * bgta
        const_9 = product_wrapper(y_p, lr, const_t, time_start)
        const_9 = product_wrapper(const_9, bgta, const_9, time_start)
        b_p = sum_wrapper(b_p, const_9, b_p, time_start)

        # b_n = b_n + lr * Y_p[i] * bgta
        const_10 = product_wrapper(y_n, lr, const_t, time_start)
        const_10 = product_wrapper(const_10, bgta, const_10, time_start)
        b_n = sum_wrapper(b_n, const_10, b_n, time_start)
            
        
        w_0_p_arr.append(w_0_p)
        w_0_n_arr.append(w_0_n)
        w_1_p_arr.append(w_1_p)
        w_1_n_arr.append(w_1_n)
        b_p_arr.append(b_p)
        b_n_arr.append(b_n)
        P_p_arr.append(P_p)
        P_n_arr.append(P_n)
        Q_p_arr.append(Q_p)
        Q_n_arr.append(Q_n)

        w_0_p_x_p_0_arr.append(w_0_p_x_p_0)
        w_1_p_x_p_1_arr.append(w_1_p_x_p_1)
        w_0_n_x_n_0_arr.append(w_0_n_x_n_0)
        w_1_n_x_n_1_arr.append(w_1_n_x_n_1)
        w_0_p_x_n_0_arr.append(w_0_p_x_n_0)
        w_0_n_x_p_0_arr.append(w_0_n_x_p_0)
        w_1_p_x_n_1_arr.append(w_1_p_x_n_1)
        w_1_n_x_p_1_arr.append(w_1_n_x_p_1)

    # finding the loss on y

    # predicting

   
    print(t, end='\r')


        



# %%
w_0 = w_0_p_arr[-1] - w_0_n_arr[-1]
w_1 = w_1_p_arr[-1] - w_1_n_arr[-1]
b = b_p_arr[-1] - b_n_arr[-1]
print(w_0, w_1, b)
