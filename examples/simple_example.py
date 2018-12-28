from cgpy.cgp import CGP, Operation

# Let's create a CGP object that represents the function
# f(x,y) = x*(y+x)-y

# Create a list of operations
operation_list = [Operation('+'), Operation('-'), Operation('*'), Operation('/')]

nr_of_variables = 2

# Each node has two (or one) inputs and one operation.
# The first nodes are the variables.
# Node 0: x
# Node 1: y
# Node 2: y+x
# Node 3: x*(y+x), which is x*node2
# Node 4: x*(y+x)-y, which is node3-y
# Note that node 0 and node 1 are automatically defined since we set the nr of variables to 2.

# The first number of each node is the operation: 0 is addition, since operation_list[0] is +.
# The second and third are the inputs, which means node 1 (y) and and node 0 (x). 
node2 = [0, 1, 0]
node3 = [2, 0, 2] # Multiplication (2) of node 0 (x) and node 2
node4 = [1, 3, 1] # Subtraction (1) of node 3 and node 1 (y).

# The last number of the gene is the index of the node that we want as output.
# In this case we want node 4 to be the output.
output_node_idx = [4]

# Add it all together and we got ourselves a gene.
gene = node2+node3+node4+output_node_idx

# That's all we need to create a CGP object.
cgp = CGP(nr_of_variables, operation_list, gene)


print("The CGP-function is:", cgp.convert2str(var_names=['x', 'y']))
print("Let x=1 and y=2. Then the function value is", cgp.eval([1,2]))
