from sympy import Matrix, Symbol, integrate, pprint

def f() -> None:
    r_0 = Matrix([600, -400, 200]) * 10**3  # km
    v_0 = Matrix([9500, 0, 0]) / 1000  # m/s -> km/s
    a_0 = Matrix([40, 0, -20]) / 1000  # m/s^2 -> km/s^2
    t = Symbol("t")

    v_t = integrate(a_0, t) + v_0
    r_t = integrate(v_t, t) + r_0

    t_value = 35 * 60 # minutes -> seconds
    v_1 = v_t.subs(t, t_value) / 10**3
    r_1 = r_t.subs(t, t_value) / 10**3

    pprint(v_t)
    pprint(r_t)
    pprint(a_0)

    pprint(f"({r_1}(10^3)) km")
    pprint(f"({v_1}(10^3)) km/s")
    pprint(f"({a_0}(10^3)) km/s^2")
    return

f()