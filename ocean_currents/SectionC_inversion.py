import matplotlib.pyplot as plt
import numpy as np
import warnings
from scipy.special import factorial, gamma, gammainc
from numpy.polynomial.polyutils import RankWarning

class Inversion:
    def __init__(self, discrete_z, discrete_U):
        self.discrete_z = discrete_z
        self.discrete_U = discrete_U

    def rms(self, x):
        return np.sqrt(np.mean(x**2))

    def standard_pedm(self, k_vect, c_til, nMax_vals=None, deltaz_T_vals=None, deltaz_B_vals=None, waterDepth=17):
        # Handling for inadequate input data.
        if len(k_vect) == 0 or len(c_til) == 0:
            return np.nan, np.nan, np.nan, np.nan, np.nan

        # Now we start with the PEDM. The steps follow those in the manuscript in section 2.1.2 Effect of limitations of measured Doppler shifts
        # Define the parameters.
        # Calculate effective depths of Doppler shift velocities based on assumption of a linear profile
        Z_eff = -(2 * k_vect) ** -1 * np.tanh(np.abs(waterDepth) * k_vect)

        # Set default values for PEDM parameter combinations if inputs left blank
        if nMax_vals is None:
            nm = min(12, round(len(k_vect) / 2))
            nMax_vals = np.arange(0, nm + 1)

        if deltaz_T_vals is None:
            depthRange = np.abs(Z_eff[0] - Z_eff[-1])
            deltaz_T_vals = np.linspace(0.01, 0.2, 20) * depthRange

        if deltaz_B_vals is None:
            depthRange = np.abs(Z_eff[0] - Z_eff[-1])
            deltaz_B_vals = np.linspace(0.02, 0.8, 20) * depthRange

        z_c = max(4 * np.min(Z_eff), -np.abs(waterDepth))  # Cutoff depth chosen as 4 times the deepest mapped depth. (Set to water depth if depth is shallower)
        # z_c is NEGATIVE by convention here

        # We loop over all PEDM parameter combinations. First initialize eps_PEDM_out as well as the other outputs.
        eps_PEDM_out = np.inf
        eps_EDM_out = np.inf

        pPEDM_out = np.nan
        pEDM_out = np.nan

        verbose = []
        combo = 0

        for nMax in nMax_vals:
            for deltaz_B in deltaz_B_vals:
                for deltaz_T in deltaz_T_vals:
                    # STEP 1: Fit the mapped Doppler shifts to a polynomial of order nMax.
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', RankWarning)
                        p1 = np.polyfit(Z_eff, c_til, nMax)

                    # STEP 2: Create additional velocity-depth pairs by linearly extrapolating up to the surface and down to cutoff depth z_c
                    zB = np.linspace(0, deltaz_B, 100)
                    zT = np.linspace(-deltaz_T, 0, 100)

                    pTop = np.polyfit(Z_eff[-1] + zT, np.polyval(p1, Z_eff[-1] + zT), 1)
                    pBottom = np.polyfit(Z_eff[0] + zB, np.polyval(p1, Z_eff[0] + zB), 1)

                    depthsExBtm = np.arange(z_c, Z_eff[0] - deltaz_B, deltaz_B)
                    depthsExTop = np.concatenate((np.arange(Z_eff[-1] + deltaz_T, 0, deltaz_T), [0]))

                    zEx = np.concatenate((depthsExBtm, Z_eff, depthsExTop))
                    cTilEx = np.concatenate(
                        (np.polyval(pBottom, depthsExBtm), c_til, np.polyval(pTop, depthsExTop)))

                    # STEP 3: Perform a second polynomial fit on the expanded set of points
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', np.RankWarning)
                        pEDM = np.polyfit(zEx, cTilEx, nMax)

                    # STEP 4: Scale polynomial coefficients defining U_EDM by n! as in equation (8) in the article.
                    pPEDM_i = pEDM / factorial(np.arange(nMax, -1, -1))

                    # STEP 5: Create a new set of linearly extrapolated points down to z_c based on the average shear of the above polynomial function in a depth interval deltaz_B/2 at the deep end of the regime.
                    zB2 = np.linspace(0, deltaz_B / 2, 100)
                    pBottom2 = np.polyfit(Z_eff[0] + zB, np.polyval(pPEDM_i, Z_eff[0] + zB), 1)

                    depthsExBtm2 = np.arange(z_c, Z_eff[0] - deltaz_B, deltaz_B / 2)

                    zEx2 = np.concatenate((depthsExBtm2, Z_eff, depthsExTop))
                    Uvals = np.concatenate((np.polyval(pBottom2, depthsExBtm2), np.polyval(pPEDM_i, Z_eff),
                                            np.polyval(pPEDM_i, depthsExTop)))

                    # STEP 6: Perform a final polynomial fit on the expanded set of points.
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', np.RankWarning)
                        pPEDM = np.polyfit(zEx2, Uvals, nMax)

                    # Calculate Doppler shifts assuming U_EDM or U_PEDM as the current profile, using the forward problem.
                    c_tilEDM = 0
                    c_tilPEDM = 0

                    if nMax < 2:
                        for n in range(nMax + 1):
                            c_tilEDM = c_tilEDM + factorial(n) * pEDM[-n - 1] * (-1 / (2 * k_vect)) ** n
                            c_tilPEDM = c_tilPEDM + factorial(n) * pPEDM[-n - 1] * (-1 / (2 * k_vect)) ** n
                    else:
                        for n in range(nMax + 1):
                            c_tilEDM = c_tilEDM + (-1 / 2) ** (n) * k_vect ** (-n) * pEDM[-n - 1] * (
                                gammainc(1 + n, -2 * k_vect * z_c)) * gamma(n + 1)
                            c_tilPEDM = c_tilPEDM + (-1 / 2) ** (n) * k_vect ** (-n) * pPEDM[-n - 1] * (
                                gammainc(1 + n, -2 * k_vect * z_c)) * gamma(n + 1)

                    # Calculate RMS differences (equation (9) in manuscript)
                    eps_EDM = self.rms(c_til - c_tilEDM)  # RMS difference for U_EDM
                    eps_PEDM = self.rms(c_til - c_tilPEDM)  # RMS difference for U_PEDM.

                    # Parameters nMax, deltaz_T, and deltaz_B are chosen to minimize eps_PEDM in practice.
                    if eps_PEDM < eps_PEDM_out:
                        eps_PEDM_out = eps_PEDM
                        pPEDM_out = pPEDM
                        nMax_out = nMax
                        deltaz_B_out = deltaz_B
                        deltaz_T_out = deltaz_T

                    if eps_EDM < eps_EDM_out:
                        eps_EDM_out = eps_EDM
                        pEDM_out = pEDM

                    combo = combo + 1

                    verbose.append({
                        'pPEDM': pPEDM,
                        'pEDM': pEDM,
                        'eps_PEDM': eps_PEDM,
                        'eps_EDM': eps_EDM,
                        'nMax': nMax,
                        'deltaz_B': deltaz_B,
                        'deltaz_T': deltaz_T
                    })
        return pPEDM_out, pEDM_out, eps_PEDM_out, eps_EDM_out, verbose

    def adjusted_pedm(self, k_vect, c_til, nMax_vals=None, deltaz_T_vals=None, deltaz_B_vals=None, waterDepth=17):
        # Added uTruth for ADCP values at discrete depths and intaking them after step 5

        # Handling for inadequate input data.
        if len(k_vect) == 0 or len(c_til) == 0:
            return np.nan, np.nan, np.nan, np.nan, np.nan

        # Now we start with the PEDM. The steps follow those in the manuscript in section 2.1.2 Effect of limitations of measured Doppler shifts
        # Define the parameters.
        # Calculate effective depths of Doppler shift velocities based on assumption of a linear profile
        Z_eff = -(2 * k_vect) ** -1 * np.tanh(np.abs(waterDepth) * k_vect)

        # Set default values for PEDM parameter combinations if inputs left blank
        if nMax_vals is None:
            nm = min(12, round(len(k_vect) / 2))
            nMax_vals = np.arange(0, nm)

        if deltaz_T_vals is None:
            depthRange = np.abs(Z_eff[0] - Z_eff[-1])
            deltaz_T_vals = np.linspace(0.01, 0.2, 20) * depthRange

        if deltaz_B_vals is None:
            depthRange = np.abs(Z_eff[0] - Z_eff[-1])
            deltaz_B_vals = np.linspace(0.02, 0.8, 20) * depthRange

        z_c = max(4 * np.min(Z_eff), -np.abs(
            waterDepth))  # Cutoff depth chosen as 4 times the deepest mapped depth. (Set to water depth if depth is shallower)
        # z_c is NEGATIVE by convention here

        # We loop over all PEDM parameter combinations. First initialize eps_PEDM_out as well as the other outputs.
        eps_PEDM_out = np.inf
        eps_EDM_out = np.inf

        pPEDM_out = np.nan
        pEDM_out = np.nan

        verbose = []
        combo = 0

        for nMax in nMax_vals:
            for deltaz_B in deltaz_B_vals:
                for deltaz_T in deltaz_T_vals:
                    # STEP 1: Fit the mapped Doppler shifts to a polynomial of order nMax.
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', np.RankWarning)
                        p1 = np.polyfit(Z_eff, c_til, nMax)

                    # STEP 2: Create additional velocity-depth pairs by linearly extrapolating up to the surface and down to cutoff depth z_c
                    zB = np.linspace(0, deltaz_B, 100)
                    zT = np.linspace(-deltaz_T, 0, 100)

                    pTop = np.polyfit(Z_eff[-1] + zT, np.polyval(p1, Z_eff[-1] + zT), 1)
                    pBottom = np.polyfit(Z_eff[0] + zB, np.polyval(p1, Z_eff[0] + zB), 1)

                    depthsExBtm = np.arange(z_c, Z_eff[0] - deltaz_B, deltaz_B)
                    depthsExTop = np.concatenate((np.arange(Z_eff[-1] + deltaz_T, 0, deltaz_T), [0]))

                    zEx = np.concatenate((depthsExBtm, Z_eff, depthsExTop))
                    cTilEx = np.concatenate((np.polyval(pBottom, depthsExBtm), c_til, np.polyval(pTop, depthsExTop)))

                    # STEP 3: Perform a second polynomial fit on the expanded set of points
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', np.RankWarning)
                        pEDM = np.polyfit(zEx, cTilEx, nMax)

                    # STEP 4: Scale polynomial coefficients defining U_EDM by n! as in equation (8) in the article.
                    pPEDM_i = pEDM / factorial(np.arange(nMax, -1, -1))

                    # STEP 5: Create a new set of linearly extrapolated points down to z_c based on the average shear of the above polynomial function in a depth interval deltaz_B/2 at the deep end of the regime.
                    zB2 = np.linspace(0, deltaz_B / 2, 100)
                    pBottom2 = np.polyfit(Z_eff[0] + zB, np.polyval(pPEDM_i, Z_eff[0] + zB), 1)

                    depthsExBtm2 = np.arange(z_c, Z_eff[0] - deltaz_B, deltaz_B / 2)

                    zEx2 = np.concatenate((depthsExBtm2, Z_eff, depthsExTop))
                    Uvals = np.concatenate((np.polyval(pBottom2, depthsExBtm2), np.polyval(pPEDM_i, Z_eff),
                                            np.polyval(pPEDM_i, depthsExTop)))

                    ''' Intakes ADCP values below a certain level
                    '''
                    indTruth = np.where(zEx2 <= np.max(self.discrete_z.flatten()))[0]
                    if len(indTruth) == 0:
                        zEx2 = np.concatenate((self.discrete_z.flatten(), zEx2))
                        Uvals = np.concatenate((self.discrete_U.flatten(), Uvals))
                    else:
                        indTruth = indTruth[-1]
                        zEx2 = np.concatenate((self.discrete_z.flatten(), zEx2[indTruth + 1:]))
                        Uvals = np.concatenate((self.discrete_U.flatten(), Uvals[indTruth + 1:]))

                    # STEP 6: Perform a final polynomial fit on the expanded set of points.
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', np.RankWarning)
                        pPEDM = np.polyfit(zEx2, Uvals, nMax)

                    # Calculate Doppler shifts assuming U_EDM or U_PEDM as the current profile, using the forward problem.
                    c_tilEDM = 0
                    c_tilPEDM = 0

                    if nMax < 2:
                        for n in range(nMax + 1):
                            c_tilEDM = c_tilEDM + factorial(n) * pEDM[-n - 1] * (-1 / (2 * k_vect)) ** n
                            c_tilPEDM = c_tilPEDM + factorial(n) * pPEDM[-n - 1] * (-1 / (2 * k_vect)) ** n
                    else:
                        for n in range(nMax + 1):
                            c_tilEDM = c_tilEDM + (-1 / 2) ** (n) * k_vect ** (-n) * pEDM[-n - 1] * (
                                gammainc(1 + n, -2 * k_vect * z_c)) * gamma(n + 1)
                            c_tilPEDM = c_tilPEDM + (-1 / 2) ** (n) * k_vect ** (-n) * pPEDM[-n - 1] * (
                                gammainc(1 + n, -2 * k_vect * z_c)) * gamma(n + 1)

                    # Calculate RMS differences (equation (9) in manuscript)
                    eps_EDM = self.rms(c_til - c_tilEDM)  # RMS difference for U_EDM
                    eps_PEDM = self.rms(c_til - c_tilPEDM)  # RMS difference for U_PEDM.

                    # Parameters nMax, deltaz_T, and deltaz_B are chosen to minimize eps_PEDM in practice.
                    if eps_PEDM < eps_PEDM_out:
                        eps_PEDM_out = eps_PEDM
                        pPEDM_out = pPEDM
                        nMax_out = nMax
                        deltaz_B_out = deltaz_B
                        deltaz_T_out = deltaz_T

                    if eps_EDM < eps_EDM_out:
                        eps_EDM_out = eps_EDM
                        pEDM_out = pEDM

                    combo = combo + 1

                    verbose.append({
                        'pPEDM': pPEDM,
                        'pEDM': pEDM,
                        'eps_PEDM': eps_PEDM,
                        'eps_EDM': eps_EDM,
                        'nMax': nMax,
                        'deltaz_B': deltaz_B,
                        'deltaz_T': deltaz_T
                    })
        return pPEDM_out, pEDM_out, eps_PEDM_out, eps_EDM_out, verbose