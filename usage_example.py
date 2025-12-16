from ocean_currents.data_loader import DataLoader
from ocean_currents.plotting import Plotter
import numpy as np

# specify a radar image sequence and a current direction
date = "20"
hour = "00"
direction = "North"

# Required Parameters
rot_angle = 4*np.pi/9 # rotation angle for nsp points

NSP_data_path = rf"C:\Users\josep\Desktop\PHYC40900_Project TP\Post_Review_Images\Generalize Process\NSP Points\NSP_{date}_{hour}.txt"
ADCP_data_path = rf"C:\Users\josep\Desktop\PHYC40900_Project TP\Post_Review_Images\Generalize Process\ADCP measurments\{date}Jan2022_{hour}00.mat"

# Load Data
loader = DataLoader(NSP_data_path, ADCP_data_path)
U_east, U_north, z = loader.load_adcp()
k_NSP, Ux_NSP, Uy_NSP = loader.load_nsp()

# # Plot
plotter = Plotter()
#plotter.plot_profiles(z, U_east)
plotter.plot_NSP(k_NSP, Ux_NSP, Uy_NSP, rot_angle)
