import sympy as sp
import numpy as np
from sympy.printing import pprint
from sympy.solvers import nsolve

'''Defining my variables.'''
t = sp.Symbol("t")
r = 2.13*t**2 - 0.0013*t**4 + 0.000034*t**4.751 # height
v = r.diff(t)

'''
r = position functon
v = velocity function
'''
r_numericalFunction = sp.lambdify(t, r, "numpy") # Numerical position function
v_numericalFunction = sp.lambdify(t, v, "numpy") # Numerical velocity function

def findApproxMaxTime():
    for i in range(0, 101, 2):
     print(f"{i}: {2.13 * i ** 2 - 0.0013 * i ** 4 + 0.000034 * i ** 4.751}") # approximately 40
    return

'''
Finding time where the rocket is at its maximum height.
'''
def calculateStuff():
    global maxTime, maxHeight
    maxTime = nsolve(v, t, 40)
    maxHeight = r.subs(t, maxTime)
    print(f"The maximum height is {maxHeight} units at {maxTime} seconds.")
    return

findApproxMaxTime()
calculateStuff()

pprint(r)
pprint(v)