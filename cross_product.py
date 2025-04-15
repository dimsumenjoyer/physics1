import numpy as np
import matplotlib.pyplot as plt

A = np.array([1, 2, 3])
B = np.array([4, 5, 6])
C = np.cross(A, B)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot vectors A, B, and C
origin = [0, 0, 0]
ax.quiver(*origin, *A, color='r', label='A', linewidth=2)
ax.quiver(*origin, *B, color='g', label='B', linewidth=2)
ax.quiver(*origin, *C, color='b', label='C', linewidth=2)

# Set plot limits
ax.set_xlim([0, max(A[0], B[0], C[0], 1)])
ax.set_ylim([0, max(A[1], B[1], C[1], 1)])
ax.set_zlim([0, max(A[2], B[2], C[2], 1)])

# Add labels and legend
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

plt.title('3D Vector Plot: A, B, and C')
plt.tight_layout()
plt.show()
