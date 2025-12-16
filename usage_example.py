from ocean_currents.data_loader import DataLoader
from ocean_currents.plotting import Plotter
from ocean_currents.SectionC_inversion import Inversion
from ocean_currents.SectionB_fitting_method import fitting_method
import numpy as np

# specify a radar image sequence and a current direction
date = "20"
hour = "00"
direction = "North"

# Required Parameters
rot_angle = 4*np.pi/9 # rotation angle for nsp points
n = 4 # polynomial degree
k_cut = 0.11
f_sigma = 1


NSP_data_path = rf"C:\Users\josep\Desktop\PHYC40900_Project TP\Post_Review_Images\Generalize Process\NSP Points\NSP_{date}_{hour}.txt"
ADCP_data_path = rf"C:\Users\josep\Desktop\PHYC40900_Project TP\Post_Review_Images\Generalize Process\ADCP measurments\{date}Jan2022_{hour}00.mat"

# Load Data
loader = DataLoader(NSP_data_path, ADCP_data_path)
U_east, U_north, z_east = loader.load_adcp()
k_NSP, Ux_NSP, Uy_NSP = loader.load_nsp()

# Section A: Plot NSP and ADCP Data
plotter = Plotter()
if direction=="East":
    plotter.plot_profiles(z_east, U_east)
else:
    plotter.plot_profiles(z_east, U_north)
Ux_transformed_NSP, Uy_transformed_NSP = plotter.plot_NSP(k_NSP, Ux_NSP, Uy_NSP, rot_angle)

# Section B: Outlier Detection and Minimisation Technique
fitting_method = fitting_method(U_east, U_north, direction, z_east)
if direction == "East":
    k_dopp_STAND, U_dopp_STAND, rmse_STAND, inter_point_STAND = fitting_method.plot_doppler_shifts(k_NSP, np.abs(Uy_transformed_NSP), k_cut, f_sigma, n)
else:
    k_dopp_STAND, U_dopp_STAND, rmse_STAND, inter_point_STAND = fitting_method.plot_doppler_shifts(k_NSP, np.abs(Ux_transformed_NSP), k_cut, f_sigma, n)

# Section C: Inversion
if direction=="East":
    discrete_U = U_east
    discrete_z = z_east
else:
    discrete_U = U_north
    discrete_z = z_east

invert = Inversion(discrete_z, discrete_U)
pPEDM_standard, pEDM_standard, a, b, c = (invert.standard_pedm(k_dopp_STAND, U_dopp_STAND))
z = np.linspace(-17, 0, 100)
U_fun_standard = lambda z: np.polyval(pPEDM_standard, z)

pPEDM_adjusted, pEDM_adjusted, d, e, f = invert.adjusted_pedm(k_dopp_STAND, U_dopp_STAND)
U_fun_adjusted = lambda z: np.polyval(pPEDM_adjusted, z)

if direction=="East":
    plotter.plot_inversion(U_fun_adjusted, U_fun_standard, z, U_east, z_east)
else:
    plotter.plot_inversion(U_fun_adjusted, U_fun_standard, z, U_north, z_east)