import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt


def main():
    filename = "test.png"

    CV_library = CV(source="RR Lyr", telescope="Kepler", nobs=10)
    lc_array = CV_library.extract_lightcurves_from_tpf()
    CV_library.plot_lightcurves(lc_array, filename)


if __name__ == "__main__":
    main()
