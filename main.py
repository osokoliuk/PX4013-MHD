import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt
import numpy as np


def main():
    filename = "test.png"
    CV_library = CV(source="KIC 4484128", telescope="Kepler", nobs=1)
    lc_array = CV_library.extract_lightcurves_from_tpf()[0].flatten(window_length=401)
    flux = np.array(lc_array.flux)
    times = (np.array(lc_array.time.unix) - np.min(lc_array.time.unix)) / 24 / 3600
    # CV_library.plot_lightcurves(lc_array, filename)
    # print(min(lc_array), max(lc_array))
    m = 64
    n = 4 * 1024
    dt = times[1] - times[0]
    t = dt * n
    print(1 / (2 * dt))
    sigmaT = 0.14
    forb = 0.5
    fourierType = "ContinuousFourier"
    spectraType = "Full"
    estimatorType = "Moment"
    s2, s3, compactS4 = CV_library.calculate_polyspectrum(
        flux, times, m, n, t, dt, sigmaT, forb, fourierType, spectraType, estimatorType
    )
    print(len(s2))
    CV_library.plot_polyspectrum(
        s2, s3, compactS4, n, t, dt, vmin=0, spectraType=spectraType
    )


if __name__ == "__main__":
    main()
