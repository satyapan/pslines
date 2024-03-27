# Horizon or Source lines plotting tool for 21-cm cosmology experiments

pslines is a Python-based tool for generating and plotting horizon lines and source lines in the cylindrical power spectra for 21-cm cosmology analyses. The tool uses the updated horizon and source line equations derived without imposing the flat sky approximation (Munshi et al. 2024, in prep) that accurately describe the signature of the horizon and the source on the power spectrum.

# Dependencies
pslines requires some standard python libraries: numpy, matplotlib, astropy, datetime.

# Installation
pslines can be imported as a python module:
```
import sys
sys.path.append('/path/to/cloned/repo/pslines')
from pslines import *
```

# Documentation
A step-by-step guide is presented in the wiki page.
