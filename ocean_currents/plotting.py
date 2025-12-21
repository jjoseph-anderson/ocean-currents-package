# plotting.py
import matplotlib.pyplot as plt
import numpy as np

class Plotter:
    def __init__(self):
        plt.rcParams.update({'font.size': 12})

    def plot_profiles(self, z, U_adcp):
        plt.plot(U_adcp, z, color='red', label='ADCP')
        plt.xlim(-0.4, 0.8)
        plt.xlabel("U [m/s]")
        plt.ylabel("Depth [m]")
        plt.legend()
        plt.grid(True)
        plt.show()

    def plot_NSP(self, k, Ux_raw, Uy_raw, rot_angle):
        Ux_transformed_NSP = np.cos(rot_angle) * Ux_raw - np.sin(rot_angle) * Uy_raw
        Uy_transformed_NSP = np.sin(rot_angle) * Ux_raw + np.cos(rot_angle) * Uy_raw

        plt.scatter(k, np.abs(Ux_transformed_NSP), label='$U_x$', color='blue')
        plt.scatter(k, np.abs(Uy_transformed_NSP), label='$U_y$', color='red')
        plt.grid(True)
        plt.legend(loc='upper right')
        plt.ylabel("U [m/s]")
        plt.xlabel("k [rad/m]")
        plt.title("Rotated NSP Doppler Shifts")
        plt.ylim(-0.05, 1.6)
        plt.xlim(0, 0.37)
        plt.show()

        return np.abs(Ux_transformed_NSP), np.abs(Uy_transformed_NSP)

    def plot_inversion(self, U_adjusted, U_standard, z, U_adcp, z_east):
        plt.plot(U_adcp, z_east, zorder=1, color='red', label='ADCP')
        plt.plot(U_adjusted(z), z, zorder=2, label="Adjusted PEDM", color='black', linestyle='--')
        plt.plot(U_standard(z), z, label="Standard PEDM", color='blue')
        plt.xlabel("U [m/s]")
        plt.ylabel("Depth [m]")
        plt.legend(loc="upper left")
        plt.grid(True)
        plt.xlim(-0.4, 0.8)
        plt.ylim(-6, 0.25)
        plt.tight_layout()
        plt.show()


