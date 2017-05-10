"""
PyPaStiX
========

"""
import ctypes
import ctypes.util

# Load the PASTIX library
libpastix_name = ctypes.util.find_library('pastix')
if libpastix_name == None:
    raise EnvironmentError("Could not find shared library: pastix."
                           "The path to libpastix.so should be in "
                           "$LIBRARY_PATH and $LD_LYBRARY_PATH.")
libpastix = ctypes.cdll.LoadLibrary(libpastix_name)

# Load the SPM library
libspm_name = ctypes.util.find_library('pastix_spm')
if libspm_name == None:
    raise EnvironmentError("Could not find shared library: spm."
                           "The path to libpastix_spm.so should be in "
                           "$LIBRARY_PATH and $LD_LYBRARY_PATH.")
libspm = ctypes.cdll.LoadLibrary(libspm_name)

__all__ = [ 'libpastix', 'libspm' ]

from .enum   import *
from .spm    import *
from .pastix import *
from .Solver import *
