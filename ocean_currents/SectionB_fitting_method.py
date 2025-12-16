import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

class fitting_method:
    def __init__(self, U_east_ADCP, U_north_ADCP, direction, z_east_ADCP):
        self.U_east_ADCP = U_east_ADCP
        self.U_north_ADCP = U_north_ADCP
        self.direction = direction
        self.z_east_ADCP = z_east_ADCP
        pass

    def dopp_shift_4poly(self, point, U_east, n):
        # Section B1 - Insert a value for the point at z=0 and then perform 4th order polynomial fit

        U_east_inter = np.insert(U_east, 0, point)
        z_east_inter = np.insert(self.z_east_ADCP, 0, 0)

        coefficients = np.polyfit(z_east_inter, U_east_inter, n)
        polynomial = np.poly1d(coefficients)

        z_values = np.linspace(-17, 0, 100)
        # values for the polynomial curve
        z_poly = np.linspace(min(z_east_inter), max(z_east_inter), 100)
        U_poly = polynomial(z_values)

        ST_east = []
        k_east = np.linspace(0.01, 0.35, 40)

        for i in range(len(k_east)):
            integral = 2 * k_east[i] * np.trapezoid(U_poly * np.e ** (2 * k_east[i] * z_poly), x=z_poly)

            # integral = ( 2 * k_east[i] )/( np.sinh(32* k_east[i]) ) * np.trapz(U_poly*np.cosh(2*k_east[i]*(16+z_poly)), x= z_poly)
            ST_east.append(integral)
        # Evaluates Stewart and Joy integral

        ST_east = np.array(ST_east)

        return k_east, ST_east, z_poly, U_poly

    def fpow(self, x, a, b, c):
        return a * x ** b + c

    def plot_doppler_shifts(self, k_data, Ux_transformed, k_cut, f_sigma, n):
        cond = (k_data > k_cut)
        filtered_k = k_data[cond]
        filtered_U = Ux_transformed[cond]

        # Perform fitting
        pars0 = (1, 0.5, 0.1)
        popt, pcov = curve_fit(self.fpow, filtered_k, filtered_U, p0=pars0, maxfev=2000)

        # Calculate residuals and plot outliers
        y_fit = self.fpow(filtered_k, *popt)
        residuals = filtered_U - y_fit
        threshold = f_sigma * np.std(residuals)

        outliers = np.abs(residuals) > threshold

        # Remove outliers and plot
        filtered_k_no_outliers = filtered_k[~outliers]
        filtered_U_no_outliers = filtered_U[~outliers]

        # Interpolate and plot
        inter_points = np.linspace(-4, 4, num=500)
        inter_k = []
        inter_U = []
        for point in inter_points:
            if self.direction == "East":
                k, U, polyk, polyU = self.dopp_shift_4poly(point, self.U_east_ADCP, n)
            else:
                k, U, polyk, polyU = self.dopp_shift_4poly(point, self.U_north_ADCP, n)
            inter_k.append(k)
            inter_U.append(U)

        # Calculate RMSE and plot the corresponding function
        rmse_values = []
        for i in range(len(inter_U)):
            f = interp1d(inter_k[i], inter_U[i], kind='linear', fill_value="extrapolate")
            interpolated_U = f(filtered_k_no_outliers)
            rmse = np.sqrt(np.mean((filtered_U_no_outliers - interpolated_U) ** 2))
            rmse_values.append(rmse)

        min_rmse_index = np.argmin(rmse_values)
        min_rmse_value = rmse_values[min_rmse_index]

        formatted_rmse = "{:.3g}".format(min_rmse_value)
        formatted_inter_point = "{:.3g}".format(inter_points[min_rmse_index])

        return inter_k[min_rmse_index], inter_U[min_rmse_index], formatted_rmse, formatted_inter_point
