import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def gaussian(x, mu, sigma):
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) / sigma)**2)

mu = 0
sigma = 1

x = np.linspace(-5, 5, 1000)

y = gaussian(x, mu, sigma)

plt.plot(x, y)
plt.title('Gaussian Function')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

result, error = quad(gaussian, -np.inf, np.inf, args=(mu, sigma))

print('Area of the Gaussian function from the bounds -5 to 5:', result)