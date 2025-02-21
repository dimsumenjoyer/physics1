from sympy import Matrix, Symbol, integrate, pprint

"""
A spaceship maneuvering near Planet Zeta is located at r⃗ =(600i^ - 400j^+200k^)(10^3)km, 
relative to the planet, and traveling at v⃗ =9500i^m/s. 
It turns on its thruster engine and accelerates with a⃗ =(40i^ - 20k^)m/s2 for 35 min.

What is the spaceship's position when the engine shuts off? Give your answer as a vector measured in km.
Express your answer in terms of the unit vectors i^,  j^, and k^.

Use the 'unit vector' button to denote unit vectors in your answer.
r⃗ = (710i^ - 400j^ + 160k^)(10^3) km
"""

def f() -> None:
    r_0 = Matrix([600, -400, 200]) * 10**3 # km
    v_0 = Matrix([9500, 0, 0]) / 1000 # m/s -> km/s
    a_0 = Matrix([40, 0, -20]) / 1000 # m/s^2 -> km/s^2
    t = Symbol("t")

    """v_t = <a_{0x}t + v_{0x}, a_{0y}t + v_{0y}, a_{0z}t + v_{0z}>"""
    v_t = integrate(a_0, t) + v_0 # symbolic
    
    """r_t = <(1/2)a_{0x}t^{2} + v_{0x}t + x_{0}, (1/2)a_{0y}t^{2} + v_{0y}t + y_{0}, (1/2)a_{0z}t^{2} + v_{0z}t + z_{0}"""
    r_t = integrate(v_t, t) + r_0 # symbolic 

    t_value = 35 * 60 # minutes -> seconds
    v_1 = (v_t.subs(t, t_value) / 10**3).evalf()
    r_1 = (r_t.subs(t, t_value) / 10**3).evalf()

    pprint(a_0); print("(10^{3}) km/s^{2}\n")
    pprint(v_t); print("(10^{3}) km/s\n")
    pprint(r_t); print("10^{3} km\n")
    
    pprint(a_0); print("(10^{3}) km/s^{2}\n")
    pprint(v_1); print("(10^{3}) km/s\n")
    pprint(r_1); print("10^{3} km\n")
    return

f()