from cgpy.cgp import CGP, Operation

from math import sqrt, pi
from random import random, seed, gauss
seed(0)

# In this case we know that the function is x-a*sin(x)-b
# but the actual values of a and b are unknown to the algorithm.

# Create a list of operations
operation_list = [Operation('+'), Operation('-'), Operation('*'), Operation('/'), Operation('sin')]

nr_of_variables = 1 # x
nr_of_parameters = 2 # a and b

# The first node is the variable x, and the second and third are the parameters a and b.

# The first number of each node is the operation: 0 is addition, since operation_list[0] is +.
# The second and third are the inputs, which means node 1 (y) and and node 0 (x). 
node3 = [4, 0, 0] # sin (4) of x (0) and then a dummy number (0) since sin only takes one input.
node4 = [2, 1, 3] # a*node3 = a*sin(x)
node5 = [1, 0, 1] # x-node4 = x - a*sin(x)
node6 = [1, 5, 2] # node5 - b = x-a*sin(x)-b

# The last number of the gene is the index of the node that we want as output.
# In this case we want node 4 to be the output.
output_node_idx = [6]

# Add it all together and we got ourselves a gene.
gene = node3+node4+node5+node6+output_node_idx

# That's all we need to create a CGP object.
cgp = CGP(nr_of_variables, operation_list, gene, nr_of_parameters=nr_of_parameters)

# Let a=0.1 and b=pi/2.
true_parameters = [0.1, pi/2.0]

# Create a list of 100 random samples of x in [0, 2*pi)
n=100
random_pnts = [[random()*2*pi] for _ in range(n)]
output_pnts = [cgp.eval(pnt, parameters=true_parameters) for pnt in random_pnts]

# And then we create a function that measures the distance between 
# the function created by the current
def error_func(cgp, random_pnts, output_pnts, parameters):
	total_error = 0.0
	for i in range(len(random_pnts)):
		pnt = random_pnts[i]
		diff = cgp.eval(pnt, parameters=parameters) - output_pnts[i]
		total_error += diff*diff
	return sqrt(total_error)

print("The CGP-function is:", cgp.convert2str(var_names=['x'], parameters=['a', 'b']))

# Create a random starting guess
parameters_guess = [random(), random()]
error_val = error_func(cgp, random_pnts, output_pnts, parameters_guess)

# Create small changes and save it if it's better.
# Rinse and repeat 10000 times.
for itr in range(10000):
	new_parameters_guess = [parameters_guess[0]+gauss(0,0.01), parameters_guess[1]+gauss(0,0.01)]
	new_error_val = error_func(cgp, random_pnts, output_pnts, new_parameters_guess)
	if new_error_val < error_val:
		error_val = new_error_val
		parameters_guess = new_parameters_guess

		print("Iteration:", itr,"  error:", error_val, " cgp function:", cgp.convert2str(var_names=['x'], parameters=parameters_guess))

