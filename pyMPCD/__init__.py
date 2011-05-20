
"""
pyMPCD - hydrodynamical simulations with MPCD
"""

# We first need to detect if we're being called as part of the numpy setup
# procedure itself in a reliable manner.
try:
    __PYMPCD_SETUP__
except NameError:
    __PYMPCD_SETUP__ = False


if __PYMPCD_SETUP__:
    import sys as _sys
    _sys.stderr.write('Running from numpy source directory.\n')
    del _sys
else:

    __all__ = ["MPCD", "MPCD_f", "test_cases"]

    from version import version as __version__

    from MPCD import MPCD_system

    import test_cases

