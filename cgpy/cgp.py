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

Let's take the the example of a 2-dimensional CGP, with
operation_table =  [+, -, *, cos]
gene = [0, 1, 0,  1, 2, 1,  3]

Call the inputs x0 and x1.
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

from random import randint, random
import sys  # Used to check python version.
from typing import List, cast  # mypy for static typing. Run 'mypy cgp.py' to get results. mypy is installable via pip.

from .operation import Operation


def _get_gene_max_values(dims: int, nr_of_parameters: int, len_of_op_table: int, nr_of_nodes: int, nodes_per_layer: int = 1) -> List[int]:
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
	assert nr_of_nodes % nodes_per_layer == 0

	dim_and_pars = dims + nr_of_parameters

	# Each node has 3 ints: the two inputs and the operation.
	len_of_gene = nr_of_nodes*3 + 1
	max_vals = [-1]*len_of_gene

	layer = 0
	for node_count in range(nr_of_nodes):
		nr_of_nodes_and_inputs_of_all_prev_layers = layer*nodes_per_layer + dim_and_pars

		max_vals[3*node_count+2] = nr_of_nodes_and_inputs_of_all_prev_layers - 1
		max_vals[3*node_count+1] = nr_of_nodes_and_inputs_of_all_prev_layers - 1
		max_vals[3*node_count] = len_of_op_table-1

		if node_count % nodes_per_layer == nodes_per_layer-1:
			layer += 1

	# The last int of the gene just points to one of the inputs, parameters or nodes and
	# calls it the outout.
	max_vals[-1] = dim_and_pars + nr_of_nodes - 1
	assert min(max_vals) >= 0
	return max_vals


def create_random_cgp(dims: int, nr_of_parameters: int, op_table: List[Operation], nr_of_nodes: int, nodes_per_layer: int = 1, fast_setup: bool = False):
	"""
	Does what it says on the tin, duh.
	"""
	len_of_op_table = len(op_table)
	max_vals = _get_gene_max_values(dims, nr_of_parameters, len_of_op_table, nr_of_nodes, nodes_per_layer=nodes_per_layer)

	random_gene = [randint(0, max_vals[i]) for i in range(len(max_vals))]

	return CGP(dims, op_table, random_gene, nr_of_parameters=nr_of_parameters, fast_setup=fast_setup)


"""
The functions _convert_cgp_2_str and _convert_rec are used by the CGP-method calc_function_str().
They create a string that is human readable and represents the mathematical function
of the CGP object.
"""
def _convert_cgp_2_str(op_table: List[Operation], gene: List[int], variable_names: List[str], nr_of_nodes: int, dims: int, parameters=[]):
	# We start at the end node and recursively work our way back.
	assert (nr_of_nodes > 0)
	assert (dims > 0)
	assert (len(variable_names) >= dims)
	assert min(gene) >= 0

	current_node_nr = gene[-1]
	return _convert_rec(op_table, gene, variable_names, nr_of_nodes, dims + len(parameters), len(parameters), current_node_nr, parameters=parameters)


def _convert_rec(op_table, gene, variable_names, nr_of_nodes, total_dims, nr_of_parameters, current_node_nr, parameters=[]):
	"""
	A recursive help function that takes a node in the CGP object and converts it into a str.
	Since a node depends on other nodes, it will recursively call this function
	but with the other nodes as input.
	If the input is a variable or a parameter, then it will simply return the
	name of the variable (as given by variable_name), or the value of the
	parameter (as given by parameters).
	"""
	assert(nr_of_parameters == len(parameters))
	# Check if the input is a variable/parameter.
	# If so, return the corresponding str.
	if(current_node_nr < total_dims):
		nr_of_vars = total_dims-nr_of_parameters
		if(current_node_nr < nr_of_vars):
			return variable_names[current_node_nr]
		else:
			return str(parameters[current_node_nr - nr_of_vars])

	op = op_table[gene[3 * (current_node_nr-total_dims) + 0]]
	nr_of_vars = total_dims-nr_of_parameters

	if op.is_binary:
		# The string is given by (str1)#(str2) if the operation is binary.
		# The #-sign is replaced by the sign of the operation. str1 and str2
		# are the strings returned by the 2 following recursive calls.
		left_str = None
		right_str = None

		# Check the left_str (str1) to see if it's a variable/parameter...
		if(gene[3 * (current_node_nr-total_dims) + 1] < total_dims):

			if gene[3 * (current_node_nr-total_dims) + 1] < nr_of_vars:
				left_str = variable_names[gene[3 * (current_node_nr-total_dims) + 1]]
			else:
				left_str = str(parameters[gene[3 * (current_node_nr-total_dims) + 1] - nr_of_vars])

		else:
			# ... if it isn't, we need a recursive function call to calc it.
			assert(gene[3 * (current_node_nr-total_dims) + 1] < current_node_nr)
			left_str = "("+_convert_rec(op_table, gene, variable_names, nr_of_nodes, total_dims, nr_of_parameters, gene[3 * (current_node_nr-total_dims) + 1], parameters=parameters)+")"

		# Repeat the process with the right_str (str2)
		if(gene[3 * (current_node_nr-total_dims) + 2] < total_dims):
			if gene[3 * (current_node_nr-total_dims) + 2] < nr_of_vars:
				right_str = variable_names[gene[3 * (current_node_nr-total_dims) + 2]]
			else:
				right_str = str(parameters[gene[3 * (current_node_nr-total_dims) + 2] - nr_of_vars])

		else:
			assert(gene[3 * (current_node_nr-total_dims) + 2] < current_node_nr)
			right_str = "("+_convert_rec(op_table, gene, variable_names, nr_of_nodes, total_dims, nr_of_parameters, gene[3 * (current_node_nr-total_dims) + 2], parameters=parameters)+")"
		return left_str+op.str+right_str
	else:  # op.is_binary

		# Things are simpler if the string is unary (non-binary). Then
		# we just return func(str1), where func is the name of the function
		# and str1 is the node that it depends on.

		middle_str = None
		# Check the middle_str (str1) to see if it's a variable/parameter...
		if(gene[3 * (current_node_nr-total_dims) + 1] < total_dims):

			if gene[3 * (current_node_nr-total_dims) + 1] < nr_of_vars:
				middle_str = variable_names[gene[3 * (current_node_nr-total_dims) + 1]]
			else:
				middle_str = str(parameters[gene[3 * (current_node_nr-total_dims) + 1] - nr_of_vars])

		else:
			# ... if it isn't, we need a recursive function call to calc it.
			middle_str = _convert_rec(op_table, gene, variable_names, nr_of_nodes, total_dims, nr_of_parameters, gene[3 * (current_node_nr-total_dims) + 1], parameters=parameters)

		# NOTE: sqr means to-the-power-of-two. This is not written
		# as sqr(x), but x^2. This is why it is treated separately.
		if op.str == "sqr":
			return "("+middle_str+")^{2}"

		return op.str+"("+middle_str+")"


class CGP():
	"""
	A Cartesian Genetic Programming object.
	This is a way of denoting a mathematical function as
	a "gene" that can be used in evolutionary optimizations.
	"""
	def __init__(self, dims: int, op_table: List[Operation], gene: List[int], nr_of_parameters: int = 0, fast_setup: bool = False, nodes_per_layer: int = 1):
		assert len(gene) > 0
		assert dims > 0
		assert len(op_table) > 0
		assert nr_of_parameters >= 0

		self.op_table = op_table  # NOTE: We only copy by reference here to speed it up a little.
		self.gene = list(gene)
		self.nr_of_parameters = nr_of_parameters
		self.dims = dims
		self.nodes_per_layer = nodes_per_layer

		self.is_constant = None
		self.has_setup_used_nodes = False

		if not fast_setup:
			self._gene_sanity_check()
			self._setup_used_nodes_list()

		self.nr_of_nodes = int((len(self.gene)-1)/3)+self.dims+self.nr_of_parameters


	def __eq__(self, other):
		# Compares all elements in one cpg to all elements in another
		# cgp to check if they are the same.
		v1 = vars(self)
		v2 = vars(other)

		assert len(v1) == len(v2)
		for key in v1:
			if v1[key] != v2[key]:
				return False
		return True


	def _gene_sanity_check(self) -> None:
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

	def _setup_used_nodes_list(self):
		# Nr of nodes excluding the output node
		nr_of_nodes = int((len(self.gene)-1)/3)+self.dims+self.nr_of_parameters
		self.has_setup_used_nodes = True

		class CGPNode():
			"""Temporary node object"""
			def __init__(self, upstream1, upstream2, is_used: int = False):
				self.upstream1 = upstream1
				self.upstream2 = upstream2
				self.is_used = is_used

			def update_is_used(self):
				self.is_used = True
				if self.upstream1 is not None:
					self.upstream1.update_is_used()
				if self.upstream2 is not None:
					assert self.upstream1 is not None
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
					assert nodes[self.gene[gene_counter]] is not None
					assert nodes[self.gene[gene_counter+1]] is not None
				else:
					nodes[i] = CGPNode(nodes[self.gene[gene_counter]], None)
					assert nodes[self.gene[gene_counter]] is not None

				gene_counter += 2
		assert gene_counter == len(self.gene)-1

		assert len(nodes) > self.gene[gene_counter]
		nodes[self.gene[gene_counter]].update_is_used()

		# See if any variables are used
		is_any_varible_used = False
		for d in range(self.dims):
			if nodes[d].is_used:
				is_any_varible_used = True
		self.is_constant = not is_any_varible_used

		# Remove the parameters and variables
		nodes = nodes[self.dims+self.nr_of_parameters:]

		assert len(nodes) == int((len(self.gene)-1)/3)
		self.used_nodes = [node.is_used for node in nodes]

	def eval(self, X: List[float], parameters: List[float] = []) -> float:
		"""
		Evaluates the function at point X using the parameters
		in parameters.
		"""
		if self.nr_of_parameters != len(parameters):
			raise Exception("The parameter dimensions are off. Expected", self.nr_of_parameters, " but got", len(parameters))

		if self.dims != len(X):
			raise Exception("There is a mismatch in dimensions in CGP. There should be", self.dims, " but the input is ", len(X))

		# Python 2 doesn't convert int/int to float, so the input
		# has to be converted to floats to make sure it is consistent.
		if sys.version_info[0] < 3:
			X = [float(val) for val in X]
			if self.nr_of_parameters:
				parameters = [float(val) for val in parameters]

		# Combined dimensionality of the variables and parameters.
		total_dims: int = len(X) + len(parameters)

		# Okay, so this is a litte weird. But n is the total number of
		# nodes used. A node is something that has a value and (possibly)
		# connections to other nodes. The returned value is the value of
		# the last node.
		n: int = int((len(self.gene)-1)/3) + total_dims
		all_node_vals: List[float] = [0.0] * n
		has_node_val_been_set: List[bool] = [False]*n

		# The inputs (variables and parameters) are the first nodes.
		for i in range(total_dims):
			if i < len(X):
				all_node_vals[i] = X[i]
				has_node_val_been_set[i] = True
			else:
				all_node_vals[i] = parameters[i-len(X)]
				has_node_val_been_set[i] = True

		# Okay, so let's step thru all the other nodes.
		node_nr = total_dims
		gene_counter = 0
		for node_nr in range(total_dims, n):
			if not self.has_setup_used_nodes or self.used_nodes[node_nr-total_dims]:
				assert(gene_counter < len(self.gene))
				assert(self.gene[gene_counter] < len(self.op_table))

				# All the nodes (except for the inputs) have an
				# operation (such as +, - cos()) and connections
				# to 1 or 2 (depending on the operation) "older"
				# nodes. "Older" means that the nodes are found
				# earlier in the list.
				op = self.op_table[self.gene[gene_counter]]
				gene_counter += 1
				node_val: float

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
				assert (not has_node_val_been_set[node_nr])
				all_node_vals[node_nr] = node_val
				has_node_val_been_set[node_nr] = True
			else:
				gene_counter += 3
		if self.has_setup_used_nodes:
			assert(sum([x is False for x in has_node_val_been_set]) == sum(x is False for x in self.used_nodes))
		else:
			assert all(has_node_val_been_set)
		assert has_node_val_been_set[self.gene[gene_counter]]
		assert gene_counter == len(self.gene)-1

		# Return the value of the last node.
		return all_node_vals[self.gene[gene_counter]]

	def calc_function_str(self, parameters=[], var_names: List[str] = None) -> str:  # TODO: Add the number of significant figures of the parameters.
		# TODO: Is this used?
		"""
		This function returns the mathematical function of the CGP as a string.
		For example 'x0+cos(0.4*x0)'.
		"""
		if len(parameters) != self.nr_of_parameters:
			raise Exception("Wrong number of parameters in the convert2str function.")

		# The variable names are automatically set to 'x1', 'x2', 'x3' and so
		# on if they are not given as inputs.
		if var_names is None:
			var_names = ["x"+str(i+1) for i in range(self.dims)]
		var_names = cast(List[str], var_names)  # Ignore this. It's only for mypy.

		if len(var_names) != self.dims:
			raise Exception("The input var_names has incorrect dimensions. Expected", self.dims, " but got", len(var_names))

		total_dims: int = len(parameters) + self.dims
		nr_of_nodes: int = int((len(self.gene)-1)/3) + total_dims

		return _convert_cgp_2_str(self.op_table, self.gene, var_names, nr_of_nodes, self.dims, parameters=parameters)

	def get_mutated_copy(self, mute_rate: float = 0.1) -> 'CGP':
		"""
		Creates a new CGP object by creating a new gene, that is a mutated
		version of the one in this object.

		# TODO: This function will be changed. In the future it will randomly
		change the gene UNTIL it alters a part of the gene that is used. A lot
		of the gene is actually not used (fun fact, the same goes for the
		human genome).
		"""
		if mute_rate < 0:
			raise Exception("The mute_rate in get_mutated_copy of CGP is less than than zero. This is not allowed.")
		elif mute_rate > 1:
			print("mute_rate in get_mutated_copy of CGP is more than 1. Setting it to 1.")
			mute_rate = 1.0

		new_gene = list(self.gene)
		nr_of_nondim_and_non_par_nodes: int = self.nr_of_nodes - self.dims - self.nr_of_parameters
		max_vals: List[int] = _get_gene_max_values(self.dims, self.nr_of_parameters, len(self.op_table), nr_of_nondim_and_non_par_nodes, nodes_per_layer=self.nodes_per_layer)

		# Mutate the gene.
		for i in range(len(new_gene)):
			if random() < mute_rate:
				new_gene[i] = randint(0, max_vals[i])

		# Let's throw in a few test cases while we are at it.
		assert len(max_vals) == len(new_gene)
		for i in range(len(new_gene)):
			assert new_gene[i] <= max_vals[i]
			assert new_gene[i] >= 0
			assert self.gene[i] <= max_vals[i]
		assert int((len(self.gene)-1)/3) == nr_of_nondim_and_non_par_nodes

		# Create the new CGP object
		new_cgp = CGP(self.dims, self.op_table, new_gene, nr_of_parameters=self.nr_of_parameters, fast_setup=not self.has_setup_used_nodes)
		return new_cgp

	def convert2str(self, parameters=[], var_names: List[str] = None):  # TODO: Add docstr
		if len(parameters) != self.nr_of_parameters:
			print("Wrong number of parameters in the convert2str function.")

		assert len(parameters) == self.nr_of_parameters

		if var_names is None:
			var_names = ["x"+str(i+1) for i in range(self.dims)]
		var_names = cast(List[str], var_names)  # Ignore this. It's only for mypy.

		total_dims = len(parameters) + self.dims
		nr_of_nodes = int((len(self.gene)-1)/3) + total_dims

		return _convert_cgp_2_str(self.op_table, self.gene, var_names, nr_of_nodes, self.dims, parameters=parameters)

	def print_function(self, parameters=[], var_names=None):  # TODO: This should be removed
		print(self.convert2str(parameters=parameters, var_names=var_names))

	def which_variables_and_parameters_are_used(self) -> List[int]:
		"""
		Returns a list of all variables and parameters that are actually
		used.
		Note that if the function is x-x, x/x, 0*x and so on, then
		the variable x IS used.

		Example: There are 3 variables (let's call them x, y, z) and two
		parameters (a & b) and the function is cos(y)*b+z/z.
		The output would be [1, 2, 4], where 1 is the y, 2 is z and 4 is b.
		"""
		# TODO: This returns a list of ints, but which_parameters_are_used returns a list of bools. Change this to bools.
		# TODO: Add some test for this one in the unit tests.

		# Combined dimensionality of the variables and parameters.
		total_dims = self.dims + self.nr_of_parameters
		n = int((len(self.gene)-1)/3) + total_dims
		node_depends_on_which_nodes: List[List[int]] = [[] for _ in range(n)]

		gene_counter = 0

		for i in range(total_dims):
			if i < self.dims:
				node_depends_on_which_nodes[i].append(i)
			else:
				node_depends_on_which_nodes[i].append(i)

		for node_nr in range(total_dims, n):
			op = self.op_table[self.gene[gene_counter]]
			gene_counter += 1

			# The node has 2 connections if the operation is binary,
			# and 1 connection otherwise.
			if op.is_binary:
				x1 = node_depends_on_which_nodes[self.gene[gene_counter]]
				gene_counter += 1
				x2 = node_depends_on_which_nodes[self.gene[gene_counter]]
				gene_counter += 1

				# Combine the two lists into one list.
				# The new list contains one instance of each value. No doubles.
				x = list(set(x2).union(set(x1)))
			else:
				x = list(node_depends_on_which_nodes[self.gene[gene_counter]])
				gene_counter += 1
				gene_counter += 1

			node_depends_on_which_nodes[node_nr] = x

		assert gene_counter == len(self.gene)-1
		return sorted(node_depends_on_which_nodes[self.gene[gene_counter]])

	def which_parameters_are_used(self) -> List[bool]:
		"""
		Returns a boolean list of length self.nr_of_parameters.
		Element i is True iff parameter nr i is used in the function.
		"""
		# TODO: This shouldn't be re-calculated every time.
		var_and_par = self.which_variables_and_parameters_are_used()

		is_parameter_used = [False for _ in range(self.nr_of_parameters)]

		for x in var_and_par:
			# Ignore the variables, and only focus on the parameters.
			if x >= self.dims:
				par_nr = x - self.dims
				is_parameter_used[par_nr] = True
		return is_parameter_used


if __name__ == '__main__':
	# This are all operations that the CGP is allowed to use. These are not set in stone.
	op_table = [Operation("+"), Operation("*"), Operation("sin"), Operation("cos"), Operation("sqr"), Operation("-"), Operation("log"), Operation("/"), Operation("sqrt")]

	dims = 2
	nr_of_parameters = 0
	nr_of_nodes = 5
	cgp = create_random_cgp(dims, nr_of_parameters, op_table, nr_of_nodes)

	print(cgp.eval([0.5, 1.5]))

	new_cgp = cgp.get_mutated_copy()
	print(new_cgp.eval([0.5, 1.5]))
