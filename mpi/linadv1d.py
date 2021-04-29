import numpy as np
import matplotlib.pyplot as plt
import os

os.system("cat ini_*.txt > ini.txt")
os.system("cat sol_*.txt> sol.txt")
ini = np.loadtxt("ini.txt")
sol = np.loadtxt("sol.txt")

plt.plot(ini[:,0], ini[:,1], label='Initial')
plt.plot(sol[:,0], sol[:,1], label='Final')
plt.legend()
plt.grid(True)
plt.show()
