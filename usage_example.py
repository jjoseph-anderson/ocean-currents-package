from ocean_currents.data_loader import DataLoader

# Load data
date = "20"
hour = "00"
direction = "North"
NSP_data_path = rf"C:\Users\josep\Desktop\PHYC40900_Project TP\Post_Review_Images\Generalize Process\NSP Points\NSP_{date}_{hour}.txt"
ADCP_data_path = rf"C:\Users\josep\Desktop\PHYC40900_Project TP\Post_Review_Images\Generalize Process\ADCP measurments\{date}Jan2022_{hour}00.mat"



loader = DataLoader(NSP_data_path, ADCP_data_path)
U_east, U_north, z = loader.load_adcp()

# Fit PEDM
# pedm = PEDM(water_depth=16)
# p = pedm.fit(k, Ux)
# U_fun = lambda z: pedm.evaluate(p, z)
#
# # Plot
# plotter = Plotter()
# plotter.plot_profiles(z, U_east, U_fun, U_fun)
