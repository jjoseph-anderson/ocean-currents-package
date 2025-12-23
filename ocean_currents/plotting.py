# plotting.py
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import pandas as pd

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

    def plot_polar_histogram(self, Du, Dm):
        fig = plt.figure(figsize=(40 / 2.54, 12 / 2.54))
        fig.subplots_adjust(wspace=0.3)
        bin_size = 4
        # Current directions
        ax1 = plt.subplot(1, 1, 1, projection='polar')
        ax1.set_theta_zero_location('N')
        ax1.set_theta_direction(-1)
        ax1.tick_params(axis='y', labelsize=8)
        plt.grid(axis='both', color='dimgrey', linestyle='-.', linewidth=0.5, alpha=0.5)

        a, b = np.histogram(Du, bins=np.arange(0, 360 + bin_size, bin_size), density=True)
        centers = np.deg2rad(np.ediff1d(b) // 2 + b[:-1])

        ax1.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, linewidth=0.7,
                color='.8', edgecolor='#ff7f0e', label='Current direction coming from')

        a, b = np.histogram(Dm, bins=np.arange(0, 360 + bin_size, bin_size), density=True)
        centers = np.deg2rad(np.ediff1d(b) // 2 + b[:-1])
        ax1.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, linewidth=0.7,
                color='.8', edgecolor='#1f77b4', label='Mean wave direction coming from ')

        plt.legend(bbox_to_anchor=(0.92, -0.1, 0., 0.), fontsize='medium', markerscale=1.2)

        ax1.set_rlabel_position(60)
        plt.tight_layout()
        plt.show()
        plt.close(fig)

    def plot_time_series(self, Du, Dm, times):
        # Plot time series
        fig, ax = plt.subplots(figsize=(10, 5))
        # Convert degrees to radians
        Du_rad = np.deg2rad(Du)
        Dm_rad = np.deg2rad(Dm)

        times1 = np.array(times)[~np.isnan((Du_rad))]
        times2 = np.array(times)[~np.isnan((Dm_rad))]

        # Remove NaN values before continuing
        Du_rad = Du_rad[~np.isnan(Du_rad)]
        Dm_rad = Dm_rad[~np.isnan(Dm_rad)]

        # Unwrap angles to remove artificial jumps
        Du_unwrapped = np.unwrap(Du_rad)
        Dm_unwrapped = np.unwrap(Dm_rad)

        # Convert back to degrees for plotting
        Du_smooth = np.rad2deg(Du_unwrapped)
        Dm_smooth = np.rad2deg(Dm_unwrapped)

        # Plot
        plt.plot(times1, np.deg2rad(Du_smooth), color='orange', label='Current direction')
        plt.plot(times2, Dm_smooth, color='blue', label='Mean wave direction')

        ax.set_ylabel('Direction [°]')
        ax.set_xlabel('Date')
        # axes[0].set_title('(a) Current Direction')
        ax.grid(True)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%dth - %H:%M'))
        plt.xticks(rotation=45)
        ax.axvline(pd.Timestamp('2022-01-18 22:00'), color='black', linestyle='-', linewidth=1)
        ax.axvline(pd.Timestamp('2022-01-20 00:00'), color='black', linestyle='-', linewidth=1)
        ax.axvspan(pd.Timestamp('2022-01-18 22:00'), pd.Timestamp('2022-01-20 00:00'), color='black', alpha=0.1,
                   label='Studied Window')
        plt.legend()
        plt.ylim(-180, 390)
        plt.tight_layout()
        plt.show()

    def plot_Hs_Tp(self, t, Tp, H):
        fig, ax1 = plt.subplots(figsize=(8, 6))

        color = 'blue'
        ax1.set_xlabel('Date')
        ax1.set_ylabel('$H_s$ [m]', color=color)
        ax1.scatter(t, H, color=color, label='$H_s$', s=12)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

        color = 'red'
        ax2.set_ylabel('$T_p$ [s]', color=color)  # we already handled the x-label with ax1
        ax2.scatter(t, Tp, color=color, label='$T_p$', s=12)
        ax2.tick_params(axis='y', labelcolor=color)

        # Set the major ticks to every 2 hours
        ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))

        # Formatting X-axis tick labels with DateFormatter
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%dth - %H:%M'))

        # Auto-rotate the date labels
        fig.autofmt_xdate()

        # Formatting
        plt.xticks(rotation=45)
        ax1.grid(True)
        plt.title("")
        # plt.legend(loc='upper left', bbox_to_anchor=(0.5, 0.5))
        ax1.axvline(pd.Timestamp('2022-01-18 22:00'), color='black', linestyle='-', linewidth=1)
        ax1.axvline(pd.Timestamp('2022-01-20 00:00'), color='black', linestyle='-', linewidth=1)
        ax1.axvspan(pd.Timestamp('2022-01-18 22:00'), pd.Timestamp('2022-01-20 00:00'), color='black', alpha=0.1,
                    label="Studied Window")
        fig.legend(loc='upper right', bbox_to_anchor=(0.88, 0.88))

        plt.show()
