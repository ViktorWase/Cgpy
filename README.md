# cgpy

A Python library for Cartesian genetic programming.

## Unit tests

Run `python -m unittest discover -v` in the root folder.

## Getting Started

Install the package using pip
`pip install cgpy`

The following code creates a random CGP and evaluates its function

```
op_table = [Operation("+"), Operation("*"), Operation("sin")]

dims = 2
nr_of_parameters = 0
nr_of_nodes = 5
cgp = create_random_cgp(dims, nr_of_parameters, op_table, nr_of_nodes)

pnt = [0.5, 1.5]
print(cgp.eval(pnt))
```

See more examples in the example folder and the documentation at www.cgpy.org

Cgpy is developed by Viktor Wase with contributions and assistance from John Brynte Turesson and Johannes Wennberg. 
