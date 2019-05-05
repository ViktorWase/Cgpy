from math import sin, cos, sqrt, log, asin, acos
from copy import copy
from typing import Union, Callable  # mypy for static typing. Run 'mypy cgp.py' to get results. mypy is installable via pip.


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
	def __init__(self, op_name: str):
		self.is_binary: bool
		self.op_name: str = op_name
		self.func: Union[Callable[[float, float], float], Callable[[float], float]]
		if op_name == "sin":
			self.func = lambda x: sin(x)
			self.is_binary = False
		elif op_name == "cos":
			self.func = lambda x: cos(x)
			self.is_binary = False
		elif op_name == "acos":
			self.func = lambda x: acos(x) if x < 1 and x > -1 else 0.0
			self.is_binary = False
		elif op_name == "asin":
			self.func = lambda x: asin(x) if x < 1 and x > -1 else 0.0
			self.is_binary = False
		elif op_name == "+":
			self.func = lambda x, y: x+y
			self.is_binary = True
		elif op_name == "-":
			self.func = lambda x, y: x-y
			self.is_binary = True
		elif op_name == "*":
			self.func = lambda x, y: x*y
			self.is_binary = True
		elif op_name == "sqr":
			self.func = lambda x: x*x
			self.is_binary = False
		elif op_name == "log":
			self.func = lambda x: log(x) if x > 0 else 0.0
			self.is_binary = False
		elif op_name == "/":
			self.func = lambda x, y: x/y if y != 0 else 0.0
			self.is_binary = True
		elif op_name == "id":
			self.func = lambda x: x
			self.is_binary = False
		elif op_name == "sqrt":
			self.func = lambda x: sqrt(x) if x > 0 else 0.0  # TODO: should this raise error instead
			self.is_binary = False
		else:
			raise ValueError('Operation not found.')
		self.str = copy(op_name)  # TODO: Is the copy actually needed?

	def __eq__(self, other):
		return self.op_name == other.op_name
