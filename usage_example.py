from ocean_currents.data_loader import DataLoader
from ocean_currents.plotting import Plotter
from ocean_currents.SectionC_inversion import Inversion
from ocean_currents.SectionB_fitting_method import fitting_method
from ocean_currents.Environmental_Conditions import read_para
from datetime import datetime
import numpy as np

# specify a radar image sequence and a current direction
#date = "20"
#hour = "00"
direction = "East"

# Required Parameters
rot_angle = 4*np.pi/9 # rotation angle for nsp points
n = 4 # polynomial degree
k_cut = 0.11
f_sigma = 1

NSP_data_path = rf"C:\Users\josep\Downloads\NSP_20_00.txt"
ADCP_data_path = rf"C:\Users\josep\Downloads\20Jan2022_0000.mat"

# Load Data
loader = DataLoader(NSP_data_path, ADCP_data_path)
U_east_init, U_north_init, z_east_init = loader.load_adcp()
U_east = U_east_init[~np.isnan(U_east_init)]
U_north = U_north_init[~np.isnan(U_north_init)]
if direction == "East":
    z_east = z_east_init[~np.isnan(U_east_init)]
else:
    z_east = z_east_init[~np.isnan(U_north_init)]

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

# Environmental Conditions
nbox = 1  # box 1, 2, 3
year = 2022
m = 1  # month
day_beg = 18
day_end = 20
dt_beg = datetime(year, m, day_beg, 0, 0)
dt_end = datetime(year, m, day_end, 2, 0)
d_beg = dt_beg.strftime("%d %b %Y")
d_end = dt_end.strftime("%d %b %Y")

finame = r'C:\\Users\\josep\\Desktop\\PHYC40900_Project TP\\2022_01_timeseries\\PAR%1i_csi_%4i%02i.txt' % (nbox, year, m)
foname = './PAR%1i_%4i%02i_directions.png' % (nbox, year, m)

lat0 = 32.07833
lon0 = 34.47233

# read WAMOS PAR_csi.txt
W_data = read_para(finame)
print(len(W_data))

Dp, Dm = [], [] # mean wave direction, peak wave direction
Du = []  # current direction
Dsea, Dswell = [], []  # wind sea and swell directions if necessary and exist

for w in W_data:
    if w.dt < dt_beg or w.dt > dt_end:
        continue
    if w.usp >= 0:
        Du.append(w.dir) # currents
    if w.dp >= 0:
        Dp.append(w.dp)  # peak dir
    if w.dm >= 0:
        Dm.append(w.dm)  # mean dir
    if w.dw >= 0:
        Dsea.append(w.dw)  # sea dir
    if w.ds >= 0:
        Dswell.append(w.ds)  # swell dir

plotter.plot_polar_histogram(Du, Dm)
