import numpy as np
import matplotlib.pyplot as plt

# Ensure that figures are not displayed on screen
plt.ioff()

# Load the data, skipping any lines that start with '#'
data = np.loadtxt('BSE[abs].out', comments='#')

# Columns: 0=E (eV), 1=Im(eps), 2=Re(eps)
E = data[:, 0]
Im_eps = data[:, 1]
Re_eps = data[:, 2]

# --- Plot E vs Im(eps) and save ---
plt.figure()
plt.plot(E, Im_eps, linewidth=1.5)
plt.xlabel('Energy (eV)')
plt.ylabel('Im(e)')
plt.title('E vs Im(e)')
plt.tight_layout()
plt.savefig('EIm.png')
plt.close()

# --- Plot E vs Re(eps) and save ---
plt.figure()
plt.plot(E, Re_eps, linewidth=1.5)
plt.xlabel('Energy (eV)')
plt.ylabel('Re(e)')
plt.title('E vs Re(e)')
plt.tight_layout()
plt.savefig('ERe.png')
plt.close()