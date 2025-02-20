import numpy as np, matplotlib.pyplot as plt

def f() -> np.ndarray:
    r_0 = np.array([600, -400, 200]) * 10**3 # km
    v_0 = np.array([9500, 0, 0]) / 1000 # m/s -> km/s
    a_0 = np.array([40, 0, -20]) / 1000 # m/s^2 -> km/s^2
    t = 35 * 60  # minutes -> seconds
    r_1 = (1/2)*a_0*t**2 + v_0*t + r_0
    return r_1 / 10**3

def plot_vectors() -> None:
    r_0 = np.array([600, -400, 200]) * 10**3  # km
    v_0 = np.array([9500, 0, 0]) / 1000  # m/s -> km/s
    a_0 = np.array([40, 0, -20]) / 1000  # m/s^2 -> km/s^2
    t = 35 * 60  # minutes -> seconds

    r_1 = (1/2)*a_0*t**2 + v_0*t + r_0
    r_1_km = r_1 / 10**3

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    ax.quiver(0, 0, 0, r_0[0]/10**3, r_0[1]/10**3, r_0[2]/10**3, color='r', label="Initial Position")
    ax.quiver(r_0[0]/10**3, r_0[1]/10**3, r_0[2]/10**3, v_0[0]*t, v_0[1]*t, v_0[2]*t, color='g', label="Velocity Vector")
    ax.quiver(r_0[0]/10**3, r_0[1]/10**3, r_0[2]/10**3, (1/2)*a_0[0]*t**2, (1/2)*a_0[1]*t**2, (1/2)*a_0[2]*t**2, color='b', label="Acceleration Effect")

    ax.scatter(r_1_km[0], r_1_km[1], r_1_km[2], color='k', marker='o', label="Final Position")

    ax.set_xlabel("X Position (km)")
    ax.set_ylabel("Y Position (km)")
    ax.set_zlabel("Z Position (km)")
    ax.set_title("Vector Representation of Motion")

    ax.legend()
    plt.show()
    return

print(f"({f()} x 10^3) km")
plot_vectors()
