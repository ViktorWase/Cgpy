"""
CGP stands for Cartesian Genetic Programming, and is a network that 
represents a function. It's really nice. Google it!

It is defined by a gene and an operation table. The operation table is 
a list of mathematical operations (for example [+, -, *, cos]). The
gene consists of n nodes and an ouput node. Each node, except for the
output, consists of 3 numbers: the operation and the two inputs. There
are always two inputs. If the operation is unary (meaning that it only
takes one input - for example cos), then the second input is simply ignored.
The last node defines which of the other nodes that is to be the output.

Let's the the example of a 2-dimensional CGP, with 
operation_table =  [+, -, *, cos]
gene = [0, 1, 0,  1, 2, 1,  3]

let's call the inputs x0 and x1.
One can illustrate what's going on by switching out the elements of the genes
that correspond to operations.
Thus [0, 1, 0,  0, 2, 1,  3] becomes [+, 1, 0,  -, 2, 1,  3]
since operation_table[0] = + and operation_table[1] = -

The other two numbers in the nodes are the indexes of the inputs. 0 in this 
case corresponds to x0, and 1 to x1. Simple. But what would 2 means? That is 
simply the first node (the one that goes: +, 1, 0,).

Thus [0, 1, 0,  0, 2, 1,  3] becomes [+, x1, x0,  -, node0, x1,  3]

The last number works on the same principle and hence 3 corresponds to
the second node.

Okay, so what does all of this mean? Well, the second node is the output.
And the second node is the first node minus x1. The first node is x1 plus x0.
Which means that [0, 1, 0,  1, 2, 1,  3] is (x1+x0)-x1.
"""

from math import sin, cos, sqrt, log, pow, exp, fabs, asin, acos, pi
from random import randint, random
from copy import copy

class Operation():
	"""
	A mathematical operation. We have a couple of the elementary ones, but 
	more might be added later. 

	Each object has:
		- a name by which it is identified
		- the function itself
		- a dual version of the function. This means 
		  that dual numbers are used.
	"""
	def __init__(self, op_name):
		self.op_name = op_name
		if op_name == "sin":
			self.func = lambda x: sin(x)
			self.is_binary = False
		elif op_name == "cos":
			self.func = lambda x: cos(x)
			self.is_binary = False
		elif op_name == "acos":
			self.func = lambda x: acos(x) if x<1 and x>-1 else 0.0
			self.is_binary = False
		elif op_name == "asin":
			self.func = lambda x: asin(x) if x<1 and x>-1 else 0.0
			self.is_binary = False
		elif op_name == "+":
			self.func = lambda x,y: x+y
			self.is_binary = True
		elif op_name == "-":
			self.func = lambda x,y: x-y
			self.is_binary = True
		elif op_name == "*":
			self.func = lambda x,y: x*y
			self.is_binary = True
		elif op_name == "sqr":
			self.func = lambda x: x*x
			self.is_binary = False
		elif op_name == "log":
			self.func = lambda x: log(x) if x>0 else 0.0
			self.is_binary = False
		elif op_name == "/":
			self.func = lambda x, y: x/y if y!=0 else 0.0
			self.is_binary = True
		elif op_name == "id":
			self.func = lambda x: x
			self.is_binary = False
		else:
			assert(False)
		self.str = copy(op_name)

def get_gene_max_values(dims, nr_of_parameters, len_of_op_table, nr_of_nodes, nodes_per_layer=1):
	"""
	A gene is a list of n ints, that define the CGP. Each such number has 
	a minimum value, and a maximum value. The minimum value is always zero.
	This function will return a list of the n maximum values.
	"""
	# Check the inputdata
	assert nodes_per_layer >= 1
	assert nr_of_nodes > 0
	assert dims > 0
	assert nr_of_parameters >= 0

	# The number of nodes has to be divisible by nodes_per_layer. 
	# Otherwise the number of layers won't be an int, and that is strange.
	assert nr_of_nodes%nodes_per_layer == 0 

	dim_and_pars = dims + nr_of_parameters

	# Each node has 3 ints: the two inputs and the operation.
	len_of_gene = nr_of_nodes*3 + 1
	max_vals = [None]*len_of_gene

	layer = 0
	for node_count in range(nr_of_nodes):
		nr_of_nodes_and_inputs_of_all_prev_layers = layer*nodes_per_layer + dim_and_pars

		max_vals[3*node_count+2] = nr_of_nodes_and_inputs_of_all_prev_layers - 1
		max_vals[3*node_count+1] = nr_of_nodes_and_inputs_of_all_prev_layers - 1
		max_vals[3*node_count] = len_of_op_table-1

		if node_count%nodes_per_layer == nodes_per_layer-1:
			layer += 1

	# The last int of the gene just points to one of the inputs, parameters or nodes and
	# calls it the outout.
	max_vals[-1] = dim_and_pars + nr_of_nodes - 1
	return max_vals

def create_random_cgp(dims, nr_of_parameters, op_table, nr_of_nodes, nodes_per_layer=1):
	"""
	Does what it says on the tin, duh.
	"""
	len_of_op_table = len(op_table)
	max_vals = get_gene_max_values(dims, nr_of_parameters, len_of_op_table, nr_of_nodes, nodes_per_layer=nodes_per_layer)

	random_gene = [randint(0, max_vals[i]) for i in range(len(max_vals))]

	return CGP(dims, op_table, random_gene, nr_of_parameters=0)

class CGP():
	"""
	A Cartesian Genetic Programming object.
	This is a way of denoting a mathematical function as
	a "gene" that can be used in evolutionary optimizations.
	"""
	def __init__(self, dims, op_table, gene, nr_of_parameters=0, fast_setup=False, nodes_per_layer=1):
		assert len(gene)>0
		assert dims > 0
		assert len(op_table) > 0
		assert nr_of_parameters >= 0
		
		self.op_table = op_table #NOTE: We only copy by reference here to speed it up a little.
		self.gene = list(gene)
		self.nr_of_parameters = nr_of_parameters
		self.dims = dims
		self.nodes_per_layer = nodes_per_layer

		self.is_constant = None

		self.has_setup_used_nodes = False

		if not fast_setup:
			self.gene_sanity_check()
			self.setup_used_nodes_list()

		self.nr_of_nodes = int((len(self.gene)-1)/3)+self.dims+self.nr_of_parameters

		assert(dims > 0)

	def gene_sanity_check(self):
		"""
		Makes sure that the gene input is consistent.
		"""
		nr_of_ins = self.dims + self.nr_of_parameters

		gene_counter = 0
		for i in range(int((len(self.gene)-1)/3)):
			assert self.gene[gene_counter] < len(self.op_table)
			gene_counter += 1
			assert self.gene[gene_counter] < i+nr_of_ins
			gene_counter += 1
			assert self.gene[gene_counter] < i+nr_of_ins
			gene_counter += 1

	def setup_used_nodes_list(self):
		# Nr of nodes excluding the output node
		nr_of_nodes = int((len(self.gene)-1)/3)+self.dims+self.nr_of_parameters
		self.has_setup_used_nodes = True
		class CGPNode():
			"""Temporary node object"""
			def __init__(self, upstream1, upstream2, is_used=False):
				self.upstream1 = upstream1
				self.upstream2 = upstream2
				self.is_used = is_used

			def update_is_used(self):
				self.is_used = True
				if self.upstream1 != None:
					self.upstream1.update_is_used()
				if self.upstream2 != None:
					assert self.upstream1 != None
					self.upstream2.update_is_used()

		nodes = [None for _ in range(nr_of_nodes)]

		gene_counter = 0
		for i in range(nr_of_nodes):
			if i < self.dims + self.nr_of_parameters:
				nodes[i] = CGPNode(None, None)
			else:
				op = self.op_table[self.gene[gene_counter]]
				gene_counter += 1

				if op.is_binary:
					nodes[i] = CGPNode(nodes[self.gene[gene_counter]], nodes[self.gene[gene_counter+1]])
					assert nodes[self.gene[gene_counter]] != None
					assert nodes[self.gene[gene_counter+1]] != None
				else:
					nodes[i] = CGPNode(nodes[self.gene[gene_counter]], None)
					assert nodes[self.gene[gene_counter]] != None

				gene_counter += 2
		assert gene_counter == len(self.gene)-1

		print(self.gene)

		assert len(nodes) > self.gene[gene_counter]
		nodes[self.gene[gene_counter]].update_is_used()

		#See if any variables are used
		is_any_varible_used = False
		for d in range(self.dims):
			if nodes[d].is_used:
				is_any_varible_used = True
		self.is_constant = not is_any_varible_used

		# Remove the parameters and variables 
		nodes = nodes[self.dims+self.nr_of_parameters:]

		assert len(nodes) == int((len(self.gene)-1)/3)
		self.used_nodes = [node.is_used for node in nodes]


	def eval(self, X, parameters = []):
		"""
		Evaluates the function at point X using the parameters
		in parameters.
		"""
		assert(self.nr_of_parameters == len(parameters))
		if self.dims != len(X):
			print("There is a mismatch in dimensions in CGP. There should be", self.dims, " but the input is ", len(X))
		assert(self.dims == len(X))

		# Combined dimensionality of the variables and parameters.
		total_dims = len(X) + len(parameters)

		# Okay, so this is a litte weird. But n is the total number of 
		# nodes used. A node is something that has a value and (possibly)
		# connections to other nodes. The returned value is the value of
		# the last node.
		n = int((len(self.gene)-1)/3) + total_dims
		all_node_vals = [None] * n 

		# The inputs (variables and parameters) are the first nodes.
		for i in range(total_dims):
			if i < len(X):
				all_node_vals[i] = X[i]
			else:
				all_node_vals[i] = parameters[i-len(X)]

		# Okay, so let's step thru all the other nodes.
		node_nr = total_dims
		gene_counter = 0
		for node_nr in range(total_dims, n):
			if self.setup_used_nodes_list==False or self.used_nodes[node_nr-total_dims]:
				assert(gene_counter<len(self.gene))
				assert(self.gene[gene_counter]<len(self.op_table))

				# All the nodes (except for the inputs) have an
				# operation (such as +, - cos()) and connections
				# to 1 or 2 (depending on the operation) "older"
				# nodes. "Older" means that the nodes are found
				# earlier in the list.
				op = self.op_table[self.gene[gene_counter]]
				gene_counter += 1
				node_val = None

				# The node has 2 connections if the operation is binary,
				# and 1 connection otherwise.
				if op.is_binary:
					x1 = all_node_vals[self.gene[gene_counter]]
					gene_counter += 1
					x2 = all_node_vals[self.gene[gene_counter]]
					gene_counter += 1

					node_val = op.func(x1, x2)
				else:
					x = all_node_vals[self.gene[gene_counter]]
					gene_counter += 1
					gene_counter += 1

					node_val = op.func(x)
				assert( all_node_vals[node_nr] == None)
				all_node_vals[node_nr] = node_val
			else:
				gene_counter += 3
		#assert(sum([x==None for x in all_node_vals])==0)
		assert(sum([x==None for x in all_node_vals]) == sum(x==False for x in self.used_nodes))
		assert all_node_vals[self.gene[gene_counter]] != None
		assert gene_counter == len(self.gene)-1

		# Return the value of the last node.
		return all_node_vals[self.gene[gene_counter]]

	def get_mutated_copy(self):
		"""
		Creates a new CGp object by creating a new gene, that is a mutated
		version of the one in this object.

		# TODO: This function will be changed. In the future it will randomly 
		change the gene UNTIL it alters a part of the gene that is used. A lot
		of the gene is actually not used (fun fact, the same goes for the 
		human genome).
		"""
		mute_rate = 0.1
		new_gene = list(self.gene)
		nr_of_nondim_and_non_par_nodes = self.nr_of_nodes - self.dims - self.nr_of_parameters
		max_vals = get_gene_max_values(self.dims, self.nr_of_parameters, len(self.op_table), nr_of_nondim_and_non_par_nodes, nodes_per_layer=self.nodes_per_layer)

		for i in range(len(new_gene)):
			if random() < mute_rate:
				new_gene[i] = randint(0, max_vals[i])

		# Let's throw in a few test cases when we are at it.
		assert len(max_vals) == len(new_gene)
		for i in range(len(new_gene)):
			assert new_gene[i] <= max_vals[i]
			assert new_gene[i] >= 0
			assert self.gene[i] <= max_vals[i]
		assert int((len(self.gene)-1)/3) ==  nr_of_nondim_and_non_par_nodes

		# Create the new CGP object
		new_cgp = CGP(self.dims, self.op_table, new_gene, nr_of_parameters=self.nr_of_parameters, fast_setup=self.setup_used_nodes_list==False)
		return new_cgp

if __name__ == '__main__':
	# This are all operations that the CGP is allowed to use. These are not set in stone.
	op_table = [Operation("+"), Operation("*"), Operation("sin"), Operation("cos"), Operation("sqr"), Operation("-"), Operation("log"), Operation("/")]

	dims = 2
	nr_of_parameters = 0
	nr_of_nodes = 5
	cgp = create_random_cgp(dims, nr_of_parameters, op_table, nr_of_nodes)

	print(cgp.eval([0.5, 1.5]))

	new_cgp = cgp.get_mutated_copy()
	print(new_cgp.eval([0.5, 1.5]))
