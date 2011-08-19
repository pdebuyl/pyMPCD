# Copyright (C) 2011 Pierre de Buyl

# This file is part of pyMPCD

# pyMPCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyMPCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyMPCD.  If not, see <http://www.gnu.org/licenses/>.

"""
pyMPCD - hydrodynamical simulations with MPCD
"""
## \namespace pyMPCD::__init__
# \brief __init__.py - This file initializes the pyMPCD module.


# We first need to detect if we're being called as part of the setup
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

