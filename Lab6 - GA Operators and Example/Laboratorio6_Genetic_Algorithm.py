import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

'''
SINGLE OBJECTIVE OPTIMIZATION GENETIC ALGORTIHM
'''


# Ogni persona della popolazione iniziale si può anche chiamare combinazione e ogni persona (combinazione) è caratterizzata
# da due parametri, x1 e x2

# Parameters
# Sia x1 che x2 hanno valore minimo 0 e max 5
xmin = np.array([0, 0])
xmax = np.array([5, 5])
n = 1000 # x1 assume 100 valori random tra 0 e 5 (lo stesso per x2)
x0 = 5 * np.random.rand(n, 2) # x0 rappresenta la popolazione iniziale (se n=1000 significa che la poplazione iniziale sarà composta da 1000 persone/combinazioni)
dimensions = x0.shape
print("Dimensions of x0:", dimensions)

par = 0
bit = 8
niteration = 1000




'''
The ga_coding function is designed to take a population matrix x, where each row represents an individual in the population, 
and each column represents a variable or dimension of the problem space. The function converts these real-valued variables into 
binary strings based on the defined number of bits (bit). The transformation scales the variables between their minimum (xmin) 
and maximum (xmax) bounds before converting them to binary.

Parameter Inputs:
x: A numpy array with shape (n_pop, n_var), where n_pop is the number of individuals in the population and n_var is the number of variables per individual.
xmin: A numpy array containing the minimum values for each variable.
xmax: A numpy array containing the maximum values for each variable.
bit: An integer indicating the number of bits used for the binary encoding of each variable.

Understanding the Matrix Dimensions:
n_pop stores the number of individuals (rows in x).
n_var stores the number of variables (columns in x).
'''
def ga_coding(x, xmin, xmax, bit):
    n_pop, n_var = x.shape
    X = np.round((x - xmin) / (xmax - xmin) * (2**bit - 1))
    B = np.array([np.binary_repr(int(xi), width=bit) for xi in X.flatten()])
    return B.reshape(n_pop, -1)



def ga_decoding(b, xmin, xmax, bit):
    n_pop = len(b)
    X = np.array([int(bi, 2) for bi in b.flatten()])
    X = X.reshape(n_pop, -1)
    x = X / (2**bit - 1) * (xmax - xmin) + xmin
    return x




def ga_crossover(xb, F, bit):
    nind = len(xb)
    ngeni = xb.shape[1]
    p = F / np.sum(F)
    ruota = np.cumsum(p)
    gen = np.random.rand(nind)
    igen = [np.searchsorted(ruota, g) for g in gen]
    cross = np.round(ngeni * np.random.rand(nind // 2)).astype(int)
    xbf = np.zeros_like(xb, dtype='<U32')  # Ensuring sufficient length for concatenation

    for i in range(nind // 2):
        p1, p2 = igen[2*i], igen[2*i+1]
        cross_point = cross[i]
        # Ensure each segment is treated as a single string
        parent1_first_part = ''.join(xb[p1, :cross_point])
        parent2_second_part = ''.join(xb[p2, cross_point:])
        parent2_first_part = ''.join(xb[p2, :cross_point])
        parent1_second_part = ''.join(xb[p1, cross_point:])

        xbf[2*i] = parent1_first_part + parent2_second_part
        xbf[2*i+1] = parent2_first_part + parent1_second_part
    
    return xbf





def ga_mutation(xb, par, bit):
    nind, ngeni = xb.shape
    pmin = min(1/nind, 1/ngeni)
    pmax = max(1/nind, 1/ngeni)
    p = par * pmin + (1 - par) * pmax
    mut = np.random.rand(nind, ngeni) < p
    for i in range(nind):
        for j in range(ngeni):
            if mut[i, j]:
                xb[i, j] = '1' if xb[i, j] == '0' else '0'
    return xb





def ga_discard(x_dis, F_dis, xson_dis, Fson_dis):
    combined = np.vstack((np.column_stack((x_dis, F_dis)), np.column_stack((xson_dis, Fson_dis))))
    sorted_indices = np.argsort(-combined[:, 2])  # Sort descending by fitness
    x_best = combined[sorted_indices][:len(x_dis)]
    return x_best[:, :2], combined[sorted_indices]





# Initialize
x_parents = x0
av_F = np.zeros(niteration)  # Store average fitness per iteration
x1_animation = np.zeros((n, niteration))  # Store X1 coordinates over iterations
x2_animation = np.zeros((n, niteration))  # Store X2 coordinates over iterations


# Genetic optimization algorithm (Single objective function)
for i in range(niteration):
    F_parents = 100 - np.sum(x_parents**2, axis=1)  # F = 0.4*GM + 0.6*ML_output (Multi objective function) Non va bene il metodo dei weights
    

    # Record positions for animation (used only for animated plotting)
    av_F[i] = np.mean(F_parents)  # Record average fitness (it's used only for plotting)
    x1_animation[:, i] = x_parents[:, 0]
    x2_animation[:, i] = x_parents[:, 1]
    
    
    bit_parents = ga_coding(x_parents, xmin, xmax, bit)
    bit_sons = ga_crossover(bit_parents, F_parents, bit)
    mutated_sons = ga_mutation(bit_sons, par, bit)
    x_sons = ga_decoding(mutated_sons, xmin, xmax, bit)
    F_sons = 100 - np.sum(x_sons**2, axis=1)
    x_best, _ = ga_discard(x_parents, F_parents, x_sons, F_sons)
    x_parents = x_best  # Optimal values (the outcome of the optimization algorithm)



# Plotting average fitness
plt.figure(figsize=(10, 5))
plt.plot(av_F, '.-')
plt.grid(True)
plt.xlabel('Iterations')
plt.ylabel('Average Fitness')
plt.title('Convergence diagram - discard BEST HALF')
plt.show()




# Animated plot
fig, ax = plt.subplots()
ax.set_xlim(0, 5)
ax.set_ylim(0, 5)
line, = ax.plot([], [], 'ro', label='Current generation')  # Initial empty plot for the current population
initial_scatter = ax.scatter(x0[:, 0], x0[:, 1], color='blue', label='Initial values')
optimal_point, = ax.plot(0, 0, 'ko', markerfacecolor='k', markersize=10, label='Optimal point (0,0)')  # Assuming the optimal point is at (0,0)
ax.legend()

def init():
    line.set_data([], [])
    return line,

def update(frame):
    line.set_data(x1_animation[:, frame], x2_animation[:, frame])
    return line,

# Adjust the interval to 500 milliseconds (0.5 seconds) per frame, enable blitting for efficient redrawing
ani = FuncAnimation(fig, update, frames=niteration, init_func=init, blit=True, interval=300)
plt.grid(True)
plt.xlabel('X1')
plt.ylabel('X2')
plt.title('Evolution of Population Over Iterations')
plt.show()





