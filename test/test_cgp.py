import unittest
from cgpy.cgp import CGP, Operation, create_random_cgp
from random import gauss, randint, random
from math import fabs, cos, sin
from copy import deepcopy

class TestCgp(unittest.TestCase):
	def setUp(self):
		pass

	def test_1(self):
		score = 1.0

		# Test 1. Evaluate CGP against a few hardcoded examples
		print("Testing CGP on hardcoded examples.")
		op_table = [Operation("+"), Operation("*"), Operation("sin"), Operation("cos"), Operation("-"), Operation("sqr")]
		genes = [[0,0,0, 0], [0,1,0, 0], [0,0,0, 0],  [0,0,0, 1],  [0,1,0, 2]]
		dim_list = [1, 2, 1, 1, 2]
		funcs = [lambda x, p: x[0], lambda x, p: x[0], lambda x, p: x[0], lambda x, p: x[0]*2, lambda x, p: x[0]+x[1]]

		n = len(dim_list)
		self.assertEquals(len(genes), n)
		self.assertEquals(len(funcs), n)
		runs = 0
		failed_runs = 0
		for i in range(n):
			for fast_setup in [False, True]:
				cgp = CGP(dim_list[i], op_table, genes[i], nr_of_parameters=0, fast_setup=fast_setup)
				for j in range(10000):
					x = [gauss(0, 100) for _ in range(dim_list[i])]
					if fabs(funcs[i](x,None)-cgp.eval(x)) > 1.0e-10:
						failed_runs += 1
					runs += 1
		score *= (1.0-(failed_runs)/(runs))
		print(100.0*(1.0-(failed_runs)/runs), "percent of the CGP cases (with parameters) worked.\n")

		# Test 2. Evaluation of CGP without parameters
		print("Testing CGP evaluation without numerical parameters")
		op_table = [Operation("+"), Operation("*"), Operation("sin"), Operation("cos"), Operation("-"), Operation("sqr")]
		failed_runs = 0.0
		runs1 = 10000
		runs2 = 10
		for i in range(runs1):
			dims = randint(2,5)
			nr_of_nodes = randint(3, 15)

			tmp_cgp = create_random_cgp(dims, 0, op_table, nr_of_nodes, nodes_per_layer=1, fast_setup=random()<0.5) # TODO: Set nodes_per_layer to a random divisor of nr_of_nodes
			#gene = create_random_gene(dims, len(op_table), nr_of_nodes)
			#tmp_cpg = CGP(dims, op_table, gene)
			x1 = None
			x2 = None
			x3 = None
			x4 = None
			x5 = None
			for j in range(runs2):
				func_str = tmp_cgp.calc_function_str()

				pnt = [ gauss(0, 10) for _ in range(dims) ]
				x1 = pnt[0]
				x2 = pnt[1]
				if dims > 2:
					x3 = pnt[2]
				if dims > 3:
					x4 = pnt[3]
				if dims > 4:
					x5 = pnt[4]

				func_str = func_str.replace('^','**')
				func_str = func_str.replace('{','(')
				func_str = func_str.replace('}',')')

				diff = fabs( tmp_cgp.eval(pnt) - eval(func_str) )

				if diff > 1.0e-10:
					failed_runs += 1.0
		
		print(100.0*(1.0-(failed_runs)/(runs1*runs2)), "percent of the CGP cases without parameters worked.\n")

		score *= (1.0-(failed_runs)/(runs1*runs2))

		# Test 3. Evaluation of CGP with parameters
		print("Testing CGP evaluation with numerical parameters")
		failed_runs = 0.0
		op_table = [Operation("+"), Operation("*"), Operation("sin"), Operation("cos"), Operation("-")]
		runs1 = 30000
		runs2 = 10
		for i in range(runs1):
			dims = randint(2,5)
			nr_of_parameters = randint(1, 3)
			all_dims = nr_of_parameters + dims

			nr_of_nodes = randint(3, 15)
			tmp_cgp = create_random_cgp(dims, nr_of_parameters, op_table, nr_of_nodes, nodes_per_layer=1, fast_setup=random()<0.5) # TODO: Set nodes_per_layer to a random divisor of nr_of_nodes
			x1 = None
			x2 = None
			x3 = None
			x4 = None
			x5 = None
			for j in range(runs2):
				pnt = [ gauss(0, 10) for _ in range(dims) ]
				params = [ gauss(0, 10) for _ in range(nr_of_parameters) ]

				func_str = tmp_cgp.calc_function_str(parameters=params)

				x1 = pnt[0]
				x2 = pnt[1]
				if dims > 2:
					x3 = pnt[2]
				if dims > 3:
					x4 = pnt[3]
				if dims > 4:
					x5 = pnt[4]

				func_str = func_str.replace('^','**')
				func_str = func_str.replace('{','(')
				func_str = func_str.replace('}',')')

				diff = fabs( tmp_cgp.eval(pnt, parameters=params) - eval(func_str) )

				if diff > 1.0e-10:
					#print("is ", tmp_cpg.eval(pnt, parameters=params), "=", func_str)
					failed_runs += 1.0

		print(100.0*(1.0-(failed_runs)/(runs1*runs2)), "percent of the CGP cases (with parameters) worked.\n")
		score *=(1.0-(failed_runs)/(runs1*runs2))

		print("Testing the mutation function.")
		run = 1000
		for _ in range(run):

			dims = randint(1, 10)
			nr_of_parameters = randint(0, 5)

			nr_of_nodes = randint(3, 20)
			npl = 1
			if random()<0.5:
				npl = randint(2, 5)
				nr_of_nodes = randint(1, 3)*npl
			cgp = create_random_cgp(dims, nr_of_parameters, op_table, nr_of_nodes, fast_setup = random()>0.5, nodes_per_layer = npl)
			cgp_old = deepcopy(cgp)
			assert cgp == cgp_old
			cgp_new = cgp.get_mutated_copy(mute_rate=random()) # This has a lot of asserts and checks in it. That's basically the unit test I guess.


		print("Total score:", round(score*10, 2),"/ 10")

if __name__ == '__main__':
    unittest.main()