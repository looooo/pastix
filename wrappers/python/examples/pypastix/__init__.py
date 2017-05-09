"""
PyPaStiX
========

"""
import ctypes
import ctypes.util

def __getnrhs(nrhs, x):
    if nrhs == -1:
        if x.ndim == 1:
            nrhs = 1
        else:
            nrhs = x.shape[1]
    return nrhs

__libpastix_name = None
__libspm_name = None
libpastix = None
libspm = None

# Load the PASTIX library
libpastix_name = ctypes.util.find_library('pastix')
if libpastix_name == None:
    raise EnvironmentError("Could not find shared library: pastix")
libpastix = ctypes.cdll.LoadLibrary(libpastix_name)

# Load the SPM library
libspm_name = ctypes.util.find_library('pastix_spm')
if libspm_name == None:
    raise EnvironmentError("Could not find shared library: spm")
libspm = ctypes.cdll.LoadLibrary(libspm_name)

__all__ = [ 'libpastix', 'libspm' ]

from .enum   import *
from .spm    import *
from .pastix import *


