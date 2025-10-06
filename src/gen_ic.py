# %%
import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# Domain parameters
# ------------------------
Nx, Ny = 1024, 512   # grid points: Nx = x-dir, Ny = y-dir
Lx, Ly = 2048, 1024  # physical dimensions

# Droplet parameters
n_rows, n_cols = 20, 10
base_radius = 24
noise_amplitude = 5

# Margin to avoid wall contact
margin = base_radius + noise_amplitude  # ensures droplet fully fits inside domain

# Initialize phase-field
phi = np.zeros((Nx, Ny), order='F')   # ensure Fortran order internally

# Generate grid coordinates (consistent with MATLAB/Fortran)
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Droplet spacing (accounting for margins)
usable_Lx = Lx - 2 * margin
usable_Ly = Ly - 2 * margin
dx = usable_Lx / n_cols
dy = usable_Ly / n_rows

# ------------------------
# Generate staggered droplets (inside margins)
# ------------------------
for i in range(n_rows):
    for j in range(n_cols):
        # Compute droplet centers within safe margins
        x_center = margin + (j + 0.5 * (i % 2)) * dx
        y_center = margin + (i + 0.5) * dy

        radius = base_radius + (np.random.rand() - 0.5) * 2 * noise_amplitude
        r = np.sqrt((X - x_center)**2 + (Y - y_center)**2)
        phi[r <= radius] = 1.0

# ------------------------
# Plot initial condition
# ------------------------
plt.contourf(X, Y, phi, levels=50, cmap='jet')
plt.colorbar(label='phi')
plt.axis('equal')
plt.title('Initial Condition: Staggered Droplets (No Wall Contact)')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()

# ------------------------
# Save in Fortran order (column-major)
# ------------------------
phi.ravel(order='F').astype('float64').tofile('phi_initial.dat')
print("Saved in Fortran order:", phi.shape)
# %%
