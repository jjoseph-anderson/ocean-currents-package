# plotting.py
import matplotlib.pyplot as plt

class Plotter:
    def __init__(self):
        plt.rcParams.update({'font.size': 12})

    def plot_profiles(self, z, U_adcp):
        plt.plot(U_adcp, z, color='red', label='ADCP')
        plt.xlabel("U [m/s]")
        plt.ylabel("Depth [m]")
        plt.legend()
        plt.grid(True)
        plt.show()
