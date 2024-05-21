import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

'''
MULTI OBJECTIVE OPTIMIZATION GENETIC ALGORTIHM

Design variables: theta1, theta2
Performance metrics: dgamma, RCH
Target values of the performance metrics: target_f1, target_f2
Objective functions: obfun1, obfun2
Fitness function: F
'''



# Data
teta1_min = -9
teta1_max = -1
teta2_min = 5
teta2_max = 15
n_points = 10 # This number represents the number of data points between the min and max of the design variables. Note that the number of all possible combination is given by n_points*n_points
teta_lbounds = np.array([teta1_min, teta2_min]) * np.pi / 180
teta_ubounds = np.array([teta1_max, teta2_max]) * np.pi / 180
theta1 = np.linspace(teta1_min, teta1_max, n_points) * np.pi / 180  # generates as many linearly spaced points between teta1_min and teta1_max as the number defined in n_points and converted from degrees to radians
theta2 = np.linspace(teta2_min, teta2_max, n_points) * np.pi / 180

DV = np.array(np.meshgrid(theta1, theta2)).T.reshape(-1, 2) # Contains n_points*n_points rows and 2 columns (theta1 and theta2)
x0 = DV   # (n_points*n_points, 2) Note that later you will see that num_individuals = n_points*n_points



# Some parameters used in the computation of the two objective functions
h1 = 0.190  # [m]
h2 = 0.270  # [m]
g = 0.050   # [m]
c = 0.870   # [m]
dz = 0.05   # [m]


# Target design variables
target_f1 = 0.100            # [m]
target_f2 = 2 * np.pi / 180  # [rad]

# Genetic Algorithm parameters
par = 0.01  # probability of mutation
bit = 8
niteration = 50  # Number of generations



x = x0
gen = 0


# The function fitness is used to assign the fitness value to each
# combination of the determined generation, this fitness is simply computed
# as the inverse of the rank number assigned to each combination after the
# computation of the Pareto optimal set.
def fitness(x, h1, h2, g, c, dz, target_f1, target_f2):  # Note that x is (n_points*n_points, 2) or we can also say (num_individuals, 2) 
    
    
    
    # ----------------------------------------------------
    # Initialization and Calculations for Each Individual
    # For each person/combination of the population, we 
    # compute all the parameters, objective functions,
    # design variables, etc.
    # ----------------------------------------------------
    
    num_individuals = x.shape[0]    # (n_points*n_points)
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
        
        
        # Objective Functions (obfun1 and obfun2) are squared error functions comparing RCH and dgamma to their 
        # respective target values (target_f1 and target_f2). The closer these values are to their targets, the 
        # lower the objective function values will be. Therefore, the objective of the optimization algorithm is 
        # to minimize the two objective functions.

        obfun1[i] = ((RCH[i] - target_f1)**2) / target_f1**2   # (num_individuals, 1)
        obfun2[i] = ((dgamma[i] - target_f2)**2) / target_f2**2  # (num_individuals, 1)
        
    # -------------------------------------------------------------------
    #  Fin of initialization and Calculations for Each Individual - FIN
    # -------------------------------------------------------------------
    
    
    
    
    
    
    # ----------------------------------------------------
    # Combine and Sort
    # ----------------------------------------------------

    # All individuals, along with their objective function values, are combined into one 
    # matrix, matr. In this matrix there are as many rows as the number of all individuals/persons/combinations, 
    # where the first column represent the first design variable (theta1), the second column represent theta2, 
    # the third column is obfun1, the fourth is obfun2 and the last column contains all zeros
    matr = np.hstack((x, obfun1, obfun2, np.zeros((num_individuals, 1)))) # (num_individuals, 5) The extra column of zeros will later hold the rank values.
    # The resulting matr matrix has the following structure:
    # [theta1, theta2, obfun1, obfun2, 0]
    # [theta1, theta2, obfun1, obfun2, 0]
    # ...


    # The matrix is sorted based on the first objective function (obfun1). This helps in the next 
    # step, where you evaluate non-dominance.
    matr_sort = matr[matr[:, 2].argsort()]  # (num_individuals, 5)
    # argsort() returns the indices that would sort this array. For example, if obfun1 values were [0.5, 0.2, 0.9], argsort() 
    # would return [1, 0, 2] because 0.2 (index 1) is the smallest, followed by 0.5 (index 0), and then 0.9 (index 2).
    # Applying argsort() on matr[:, 2] gives us the indices that sort obfun1 in ascending order. Using these indices to index 
    # matr (matr[matr[:, 2].argsort()]) sorts the entire matr matrix based on the obfun1 values. This sorted matrix, matr_sort, 
    # is structured such that rows are ordered by increasing obfun1.
    # Why Sorting by obfun1?
    # Sorting the matrix by obfun1 is a preparatory step for the non-dominated sorting that follows. Here's why it's useful:
    # 1) Easier Identification of Non-Dominated Individuals: Starting with the smallest obfun1 values simplifies the process 
    # of finding non-dominated individuals. When you have the smallest obfun1 at the top, you can more easily check other 
    # objective values (obfun2 in this case) to determine non-dominance.
    # 2) Efficient Ranking: By sorting, you create a more structured approach to assigning ranks. You can start with the 
    # best (lowest) values of obfun1 and then determine the rank based on obfun2.
    # Assume you have three individuals with the following values:
    # Individual   theta1	theta2	obfun1	obfun2
    #    1	         0.3	 0.6	 0.9	 0.4
    #    2	         0.1	 0.2	 0.2	 0.5
    #    3	         0.5	 0.4	 0.5	 0.3
    
    # Matrix Before Sorting (matr):
    # [[0.3, 0.6, 0.9, 0.4, 0],
    #  [0.1, 0.2, 0.2, 0.5, 0],
    #  [0.5, 0.4, 0.5, 0.3, 0]]
    # As you can see the sequence of obfun1 is 0.9, 0.2, 0.5
    
    # Sorted Matrix (matr_sort) based on obfun1:
    # [[0.1, 0.2, 0.2, 0.5, 0],
    #  [0.5, 0.4, 0.5, 0.3, 0],
    #  [0.3, 0.6, 0.9, 0.4, 0]]
    # Now after the sorting, the sequence of obfun1 is 0.2, 0.5, 0.9
    
    # -------------------------------------------------------------------
    #  Fin of Combine and Sort - FIN
    # -------------------------------------------------------------------
    
    
    
    
    
    
    # ----------------------------------------------------
    # Non-Dominated Sorting
    # ----------------------------------------------------
    
    # This block implements non-dominated sorting. An individual is considered non-dominated if no other individual in the 
    # population is better in all objectives.
    matr_rank = []  # Pre allocate memory to the matr_rank matrix
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

        # Sorting by Ranks: Individuals that are non-dominated get the lowest rank (best rank), and are then removed from the 
        # sorting pool. This process repeats until all individuals are ranked.
        rank = np.full((len(index), 1), z + 1)
        matr_rank.append(np.hstack((matr_sort[index, :4], rank)))
        matr_sort = np.delete(matr_sort, index, axis=0)
        z += 1
        
    # -------------------------------------------------------------------
    #  Fin of Non-Dominated Sorting - FIN
    # -------------------------------------------------------------------





    # -------------------------------------------------------------------
    # Calculate Fitness
    # -------------------------------------------------------------------
    
    # Final Fitness: Here, the fitness is inversely proportional to the rank. A lower rank (more non-dominance) results 
    # in a higher fitness.
    if matr_rank:
        matr_rank = np.vstack(matr_rank)

    F = 1 / matr_rank[:, 4]  # Compute fitness as the inverse of the rank
    
    # -------------------------------------------------------------------
    # Fin of Calculate Fitness - FIN
    # -------------------------------------------------------------------
    
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




def ga_mutation(xb, par, bit):
    nind, ngeni = xb.shape
    pmin = min(1/nind, 1/ngeni)
    pmax = max(1/nind, 1/ngeni)
    p = par * pmin + (1 - par) * pmax  # Mutation probability calculation

    # Generate a random array the same shape as xb
    mut = np.random.rand(nind, ngeni)
    imut = mut < p  # Determine which genes will mutate

    # Perform mutation: flip bits where mutation occurs
    xb[imut] = 1 - xb[imut]  # Flipping bits: if 1, make it 0; if 0, make it 1

    return xb




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



# Genetic optimization algorithm
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



# # Plotting results
# plt.figure()
# plt.plot(x[:, 0], x[:, 1], 'ro', label='Pareto-Optimal points')
# plt.plot(x0[:, 0], x0[:, 1], 'bo', label='Initial values')
# plt.grid(True)
# plt.xlabel('X_1')
# plt.ylabel('X_2')
# plt.title('Design Variable Space')
# plt.legend()

# plt.figure()
# plt.plot(f1initial, f2initial, 'bo', label='Initial values')
# plt.plot(f1parent, f2parent, 'ro', label='Pareto-optimal points derived from GA')
# plt.grid(True)
# plt.xlabel('f1')
# plt.ylabel('f2')
# plt.title('Objective Functions Space')
# plt.legend()

# plt.figure()
# plt.bar(range(gen), NIR1 / len(F), color='b')
# plt.ylabel('Number of individuals with rank 1 over number of total individuals')
# plt.xlabel('Generations')
# plt.grid(True)
# plt.title('Rank behaviour')
# plt.show()







