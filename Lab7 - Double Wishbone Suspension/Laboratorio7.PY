import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

'''
MULTI OBJECTIVE OPTIMIZATION GENETIC ALGORTIHM
'''



# Data
teta1_min = -9
teta1_max = -1
teta2_min = 5
teta2_max = 15
teta_lbounds = np.array([teta1_min, teta2_min]) * np.pi / 180
teta_ubounds = np.array([teta1_max, teta2_max]) * np.pi / 180
theta1 = np.linspace(teta1_min, teta1_max, 20) * np.pi / 180
theta2 = np.linspace(teta2_min, teta2_max, 20) * np.pi / 180

DV = np.array(np.meshgrid(theta1, theta2)).T.reshape(-1, 2)
x0 = DV

h1 = 0.190  # [m]
h2 = 0.270  # [m]
g = 0.050   # [m]
c = 0.870   # [m]
dz = 0.05   # [m]

target_f1 = 0.100       # [m]
target_f2 = 2 * np.pi / 180  # [rad]

# Genetic Algorithm parameters
par = 0.03  # probability of mutation
bit = 8

niteration = 200

x = x0
gen = 0

# Example functions placeholders (you should implement these based on your MATLAB code logic)
def fitness(x, h1, h2, g, c, dz, target_f1, target_f2):
    num_individuals = x.shape[0]
    d = np.zeros((num_individuals, 1))
    h = np.zeros((num_individuals, 1))
    RCH = np.zeros((num_individuals, 1))
    dgamma = np.zeros((num_individuals, 1))
    obfun1 = np.zeros((num_individuals, 1))
    obfun2 = np.zeros((num_individuals, 1))

    for i in range(num_individuals):
        theta1 = x[i, 0]
        theta2 = x[i, 1]

        d[i] = h2 / (np.tan(theta1 + theta2) - np.tan(theta1))
        h[i] = h1 - d[i] * np.tan(theta1)

        RCH[i] = c * h[i] / (d[i] + g)
        dgamma[i] = dz / (d[i] + g)

        obfun1[i] = ((RCH[i] - target_f1)**2) / target_f1**2
        obfun2[i] = ((dgamma[i] - target_f2)**2) / target_f2**2

    matr = np.hstack((x, obfun1, obfun2, np.zeros((num_individuals, 1))))
    matr_sort = matr[matr[:, 2].argsort()]  # Sort by obfun1
    
    matr_rank = []
    z = 0
    while matr_sort.shape[0] > 0:
        index = []
        for i in range(matr_sort.shape[0]):
            dominated = False
            for j in range(matr_sort.shape[0]):
                if matr_sort[j, 3] < matr_sort[i, 3] and matr_sort[j, 2] < matr_sort[i, 2]:
                    dominated = True
                    break
            if not dominated:
                index.append(i)

        rank = np.full((len(index), 1), z + 1)
        matr_rank.append(np.hstack((matr_sort[index, :4], rank)))
        matr_sort = np.delete(matr_sort, index, axis=0)
        z += 1

    if matr_rank:
        matr_rank = np.vstack(matr_rank)

    F = 1 / matr_rank[:, 4]  # Compute fitness as the inverse of the rank
    return F, matr_rank, obfun1.flatten(), obfun2.flatten()


def ga_coding(x, xmin, xmax, bit):
    n_pop, n_var = x.shape
    b = np.zeros((n_pop, n_var * bit), dtype=int)

    for i in range(n_var):
        # Normalize and scale each variable to binary
        norm = (x[:, i] - xmin[i]) / (xmax[i] - xmin[i])
        scaled = np.floor(norm * (2 ** bit - 1)).astype(int)

        # Convert scaled integers to binary
        for j in range(bit):
            b[:, i * bit + j] = (scaled >> j) & 1

    return b



def ga_crossover(b, F, bit):
    n_pop, n_genes = b.shape
    # Ensure that crossover happens only if there is a population to work with
    if n_pop < 2:
        return b

    crossed_b = np.empty_like(b)
    for i in range(0, n_pop, 2):
        if i + 1 >= n_pop:  # If odd number of individuals, just copy the last one without crossover
            crossed_b[i, :] = b[i, :]
            break
        
        # Random crossover point
        point = np.random.randint(0, n_genes)
        
        # Perform crossover
        crossed_b[i, :point] = b[i, :point]
        crossed_b[i, point:] = b[i + 1, point:]
        crossed_b[i + 1, :point] = b[i + 1, :point]
        crossed_b[i + 1, point:] = b[i, point:]

    return crossed_b


def ga_mutation(xbf_cross, par, bit):
    # Implement your mutation function here
    pass

def ga_decoding(b, xmin, xmax, bit):
    n_pop = b.shape[0]
    n_var = len(xmin)
    x = np.zeros((n_pop, n_var))
    
    for i in range(n_var):
        # Convert binary strings back to decimal and then scale to original range
        x[:, i] = ((b[:, i*bit:(i+1)*bit] @ (2**np.arange(bit)[::-1])) / (2**bit - 1)) * (xmax[i] - xmin[i]) + xmin[i]

    return x


def ga_discard(xparent, Fparent, xson, Fson):
    # Combine parents and sons
    combined_x = np.vstack((xparent, xson))
    combined_F = np.hstack((Fparent, Fson))

    # Create a structured array to keep x values paired with their fitness
    dtype = [('x', float, (2,)), ('F', float)]
    combined = np.empty(combined_x.shape[0], dtype=dtype)
    combined['x'] = combined_x
    combined['F'] = combined_F

    # Sort the combined array by fitness in descending order
    combined_sorted = np.sort(combined, order='F')[::-1]

    # Select the top half
    midpoint = len(combined_sorted) // 2
    x_best = combined_sorted['x'][:midpoint]
    pop_ord = combined_sorted['x']

    return x_best, pop_ord




for ii in range(niteration):
    F, matr_rank, obfun1, obfun2 = fitness(x, h1, h2, g, c, dz, target_f1, target_f2)
    gen += 1

    if gen == 1:
        xinitial = DV
        f1initial = obfun1
        f2initial = obfun2

    xparent = matr_rank[:, :2]
    f1parent = matr_rank[:, 2]
    f2parent = matr_rank[:, 3]
    check = np.where(F == 1)[0]
    NIR1 = len(check)

    if np.all(F == 1):
        break

    b = ga_coding(xparent, teta_lbounds, teta_ubounds, bit)
    xbf_cross = ga_crossover(b, F, bit)
    xm = ga_mutation(xbf_cross, par, bit)
    xson = ga_decoding(xbf_cross, teta_lbounds, teta_ubounds, bit)
    Fson, matr_rank_son, _, _ = fitness(xson, h1, h2, g, c, dz, target_f1, target_f2)
    x_best, pop_ord = ga_discard(xparent, F, xson, Fson)
    x = x_best

# Plotting results
plt.figure()
plt.plot(x[:, 0], x[:, 1], 'ro', label='Pareto-Optimal points')
plt.plot(x0[:, 0], x0[:, 1], 'bo', label='Initial values')
plt.grid(True)
plt.xlabel('X_1')
plt.ylabel('X_2')
plt.title('Design Variable Space')
plt.legend()

plt.figure()
plt.plot(f1initial, f2initial, 'bo', label='Initial values')
plt.plot(f1parent, f2parent, 'ro', label='Pareto-optimal points derived from GA')
plt.grid(True)
plt.xlabel('f1')
plt.ylabel('f2')
plt.title('Objective Functions Space')
plt.legend()

plt.figure()
plt.bar(range(gen), NIR1 / len(F), color='b')
plt.ylabel('Number of individuals with rank 1 over number of total individuals')
plt.xlabel('Generations')
plt.grid(True)
plt.title('Rank behaviour')
plt.show()







