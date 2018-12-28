from cgpy.cgp import CGP, Operation, create_random_cgp
from math import sqrt
from random import gauss, seed

seed(0)

# Let's do some basic machine learning
# We create 100 random points in 2D space
n=100
random_points = [[gauss(0, 1), gauss(0,1)] for _ in range(n)]

# Then we put them through a function. Let's say x+y*y
f = lambda pnt: pnt[0]+pnt[1]*pnt[1]
out_pnts = [f(pnt) for pnt in random_points]

# And let's define a function that measures the differences
# between the cgp and the inp_pnts and out_pnts.
def error_function(cgp, inp_pnts, out_pnts):
	total_error = 0.0
	for i in range(len(inp_pnts)):
		pnt = inp_pnts[i]
		diff = cgp.eval(pnt)-out_pnts[i]

		total_error += diff*diff
	return sqrt(total_error)

# Create a list of operations
operation_list = [Operation('+'), Operation('-'), Operation('*'), Operation('/')]

# We use two variables and no parameters. Check out parameter_example.py to learn about them.
nr_of_variables = 2
nr_of_parameters = 0
len_of_op_table = len(operation_list)

# The number of nodes is up to the user. Shorter means simpler functions, and longer means more complex.
nr_of_nodes = 5

# Create an initial random cgp object.
cgp = create_random_cgp(nr_of_variables, nr_of_parameters, operation_list, nr_of_nodes)
optimization_value = error_function(cgp, random_points, out_pnts)

for itr in range(10000):
	# Mutate it and keep it if it's better than the old one.
	new_cgp = cgp.get_mutated_copy()
	new_opt_value = error_function(new_cgp, random_points, out_pnts)

	if new_opt_value<optimization_value:
		optimization_value = new_opt_value
		cgp = new_cgp
		print("Iteration:", itr,"  error:", optimization_value, " cgp function:", cgp.convert2str(var_names=['x', 'y']))