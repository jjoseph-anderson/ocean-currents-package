## Ocean-currents Package 
___
**Author:** Joseph Anderson 

This Python package has been created as part of our recent paper *Advances in extracting current profiles from X-band radar images with a focus on retrieving subsurface current* which is now available on google scholar. It's purpose is to:
- Find the Doppler shifts that are dependent on the wavenumber using the Normalized Scalar Product (NSP) method
- Define a Doppler shift curve by:
  - Performing outlier detection
  - Utilizing a minimization technique to fit to the ADCP profile
- Invert the Doppler shift curve to produce the depth-dependent current profile

### Import required packages
___

'''python 
import numpy as np
import matplotlib.pyplot as plt
import h5py
