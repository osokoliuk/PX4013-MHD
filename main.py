import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt
import numpy as np


def main():
    filename = "test.png"
    CV_library = CV(source="KIC 5520878", telescope="Kepler", nobs=1)
    lc_array = CV_library.extract_lightcurves_from_tpf()[0]
    flux = np.array(lc_array.flux)
    times = np.array(lc_array.time.unix)
    # CV_library.plot_lightcurves(lc_array, filename)
    # print(min(lc_array), max(lc_array))
    # x = np.linspace(0.1, 1e4 - 0.1, 100000)
    # lc_array = np.sin(x) + np.cos(x) ** 2 + np.sin(x) ** 3
    m = 10
    n = 256
    dt = times[0] - times[1]
    print(times)
    t = dt * n
    sigmaT = 0.14
    forb = 0.5
    spectraType = "Full"
    estimatorType = "Moment"
    s2, s3, compactS4 = CV_library.calculate_polyspectrum(
        flux, times, m, n, t, dt, sigmaT, forb, spectraType, estimatorType
    )
    CV_library.plot_polyspectrum(
        s2, s3, compactS4, n, t, vmin=0, spectraType=spectraType
    )


if __name__ == "__main__":
    main()
