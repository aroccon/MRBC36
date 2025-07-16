# %%
import numpy as np
import matplotlib.pyplot as plt

nx, ny = 128, 64

data = np.fromfile('output/p_00010000.dat', dtype=np.float64)  # or float64
print("Data size:", data.size)
data = data.reshape((ny, nx))  # rows = y, columns = x

plt.imshow(data, cmap='jet', origin='lower', aspect='auto')
plt.colorbar(label='Variable')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# %%
