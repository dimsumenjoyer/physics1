from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, mu, sigma):
    """
    Parameters:
    x: float, input value.
    mu: float, mean of the distribution.
    sigma: float, standard deviation of the distribution.
    """
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) / sigma)**2)

def riemann_sum(gaussian_func, a, b, mu, sigma, n):
    """
    Parameters:
    gaussian_func: function, the Gaussian function to integrate.
    a: float, the lower bound of the integration.
    b: float, the upper bound of the integration.
    mu: float, the mean of the Gaussian distribution.
    sigma: float, the standard deviation of the Gaussian distribution.
    n: int, the number of sub-intervals to use in the approximation.
    """
    dx = (b - a) / n  # Width of each sub-interval
    area = 0
    for i in range(n):
        x = a + i * dx  # Left endpoint of each sub-interval
        area += gaussian_func(x, mu, sigma) * dx
    return area

def trapezoidal_approx(gaussian_func, a, b, mu, sigma, n):
    """
    Parameters:
    gaussian_func: function, the Gaussian function to integrate.
    a: float, the lower bound of the integration.
    b: float, the upper bound of the integration.
    mu: float, the mean of the Gaussian distribution.
    sigma: float, the standard deviation of the Gaussian distribution.
    n: int, the number of sub-intervals to use in the approximation.
    """
    x = np.linspace(a, b, n)
    y = gaussian_func(x, mu, sigma)
    area = np.trapz(y, x)
    return area

def definite_integral(gaussian_func, a, b, mu, sigma):
    """
    Parameters:
    gaussian_func: function, the Gaussian function to integrate.
    a: float, the lower bound of the integration.
    b: float, the upper bound of the integration.
    mu: float, the mean of the Gaussian distribution.
    sigma: float, the standard deviation of the Gaussian distribution.
    """
    return quad(gaussian_func, a, b, args = (mu, sigma))[0]  # Return only the area, omit error

def plot_gaussian(gaussian_func, mu, sigma, a, b):
    """
    Parameters:
    gaussian_func: function, the Gaussian function to plot.
    mu: float, the mean of the Gaussian distribution.
    sigma: float, the standard deviation of the Gaussian distribution.
    a: float, the lower bound of the area to shade.
    b: float, the upper bound of the area to shade.
    """
    x = np.linspace(-5, 5, 1000)
    y = gaussian_func(x, mu, sigma)

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='Gaussian Function', color='blue')

    # Fill the area under the curve from a to b
    plt.fill_between(x, y, where = (x >= a) & (x <= b), color = 'skyblue', alpha = 0.5, label = f'Area from {a} to {b}')

    plt.title('Gaussian Function and Area Estimations')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axhline(0, color = 'black', lw = 0.5, ls = '--')
    plt.axvline(0, color = 'black', lw = 0.5, ls = '--')
    plt.legend()
    plt.grid()
    plt.show()

riemann_area = riemann_sum(gaussian, -5, 5, 0, 1, 1000)
trapezoidal_area = trapezoidal_approx(gaussian, -5, 5, 0, 1, 1000)
exact_area = definite_integral(gaussian, -5, 5, 0, 1)

print(f"Estimated area using Riemann sum: {riemann_area}")
print(f"Estimated area using trapezoidal approximation: {trapezoidal_area}")
print(f"Exact area using integration: {exact_area}")

riemann_sum_margin_of_error = (abs(riemann_area - exact_area) / exact_area) * 100
trapezoidal_approx_margin_of_error = (abs(trapezoidal_area - exact_area) / exact_area) * 100

print(f"Margin of error for Riemann sum: {riemann_sum_margin_of_error: 0.50f}%")
print(f"Margin of error for trapezoidal approximation: {trapezoidal_approx_margin_of_error: 0.50f}%")

plot_gaussian(gaussian, 0, 1, -5, 5)