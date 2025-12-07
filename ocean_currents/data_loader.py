# data_loader.py
import numpy as np
import h5py
from scipy.io import loadmat

class DataLoader:
    def __init__(self, nsp_path, adcp_path):
        self.nsp_path = nsp_path
        self.adcp_path = adcp_path

    def load_nsp(self):
        data = np.loadtxt(self.nsp_path, delimiter=',', skiprows=1)
        k, Ux, Uy = data[:,0], data[:,1], data[:,2]
        return k, Ux, Uy

    def load_adcp(self):
        try:
            with h5py.File(self.adcp_path, 'r') as hdf:
                U_east = np.array(hdf['East_vel']).flatten()
                U_north = np.array(hdf['North_vel']).flatten()
                z = np.array(hdf['Param/z_vec']).flatten()
        except OSError:
            data = loadmat(self.adcp_path)
            U_east = data['East_vel'][0]
            U_north = data['North_vel'][0]
            z = data['Param']['z_vec'].flatten()[0][0]
        return U_east, U_north, z
