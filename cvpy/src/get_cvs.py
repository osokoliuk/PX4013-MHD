"""
Module      : cvpy.get_cvs.py
Description : Parser for cataclysmic variable lightcurves
Copyright   : (c) Oleksii Sokoliuk, 2026
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

This script finds and downloads all of the available lightcurves from Kepler,
TESS and K2 missions for a given source coordinates.
"""

import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk

ax = plt.subplot(111)

source = "V344 Lyr"

search_result = lk.search_targetpixelfile(source, author="Kepler")
search_download = search_result.download_all()[:5]
for dataset in search_download:
    lc = dataset.to_lightcurve()
    lc.plot(c="k", ax=ax)
plt.savefig("test.png")
