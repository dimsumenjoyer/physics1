import typing

def trapz(f: typing.Callable[[float], float], a: int, b: float, dx: float) -> float:
    n: int = int((b - a) / dx)
    dx: float = (b - a) / n
    integral: float = f(a) + f(b)

    for k in range(1, n):
        x_k: float = a + k*dx
        integral += 2 * f(x_k)

    return (dx / 2) * integral

print(trapz(lambda x: x**2, 0, 2, 0.5))
