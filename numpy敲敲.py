import numpy as np
from matplotlib import pyplot as plt

x = np.arange(-11, 11)
y = 2 * x**2 + 5*x + 5
z = 3 * y
plt.title("Matplotlib demo")
plt.xlabel("x axis caption")
plt.ylabel("y axis caption")
plt.zlabel("z axis caption")
plt.plot(x, y, z)
plt.show()
'''
a = np.array([1, 2, 3], dtype=complex, copy=True, ndmin=2)
print(a)

b = np.arange(24)
print(b)
print(b.ndim)
print(a.ndim)
c = b.reshape(2, 4, 3)
print(c.ndim)
print(c)
'''