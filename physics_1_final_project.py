import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy.printing import pprint
from sympy.solvers import nsolve

t = sp.Symbol("t")
r = 2.13*t**2 - 0.0013*t**4 + 0.000034*t**4.751 # height
v = r.diff(t)

pprint(r)
pprint(v)

r_numericalFunction = sp.lambdify(t, r, "numpy")  # Numerical position function
v_numericalFunction = sp.lambdify(t, v, "numpy")  # Numerical velocity function

time_values = np.arange(0, 101, 2) # Generate time values (0 to 100 seconds in 2 second intervals)
position_values = r_numericalFunction(time_values) 
velocity_values = v_numericalFunction(time_values)

def findApproxMaxTime():
    for i in range(0, 101, 2):
     print(f"{i}: {2.13 * i ** 2 - 0.0013 * i ** 4 + 0.000034 * i ** 4.751}") # approximately 40
    return

def calculateStuff():
    global maxTime, maxHeight
    maxTime = nsolve(v, t, 40)
    maxHeight = r.subs(t, maxTime)
    print(f"The maximum height is {maxHeight} units at {maxTime} seconds.")
    return

def plot():
    plt.figure(figsize = (10, 6))

    plt.plot(time_values, position_values, label = "Position (Height)", color = "blue")
    plt.axvline(x = float(maxTime), color = "blue", linestyle = "--", label = f"Max Height at t = {float(maxTime):}")

    plt.plot(time_values, velocity_values, label = "Velocity", color = "red")
    plt.axhline(y = 0, color = "black", linestyle = "--", label = "Velocity = 0")

    plt.title("Position & Velocity of Rocket vs. Time")
    plt.xlabel("Time (Seconds)")
    plt.ylabel("Value")
    plt.legend()
    plt.grid(True)

    plt.show()
    return

findApproxMaxTime()
calculateStuff()
plot()