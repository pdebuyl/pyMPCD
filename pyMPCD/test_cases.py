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
test_cases.py - This file contains test cases for the pyMPCD simulation package.
"""
## \namespace pyMPCD::test_cases
# \brief test_cases.py - This file contains test cases for the pyMPCD simulation package.
#
# These test_cases are intented to demonstrate some of the capabilities of the 
# package and to provide initialization of tests.
#
# \b Example of usage
# \code
#import pyMPCD
#box = pyMPCD.test_cases.eight_PBC()
#for i in range(10):
#    box.one_full_step()
# \endcode

from pyMPCD import MPCD_system
import numpy as np

## Returns a 8**3 system with density 10 and temperature T=.33.
# Periodic boundary conditions in x,y and z.
def eight_PBC():
    """
    Returns a 8**3 system with density 10 and temperature T=.33.
    Periodic boundary conditions in x,y and z.
    """
    a = MPCD_system( (8,8,8) , 10, 1.)
    a.null_shift()
    print a
    a.tau = .5
    a.init_r()
    a.init_v(.33)
    return a

## Returns a 8**3 system with density 10 and temperature T=.33.
# Periodic boundary conditions in y and z.
# Thermostatted walls in direction x.
def eight_PBC_wallx():
    """
    Returns a 8**3 system with density 10 and temperature T=.33.
    Periodic boundary conditions in y and z.
    Thermostatted walls in direction x.
    """
    a = MPCD_system( (8,8,8) , 10, 1.)
    a.null_shift()
    print a
    a.tau = .5
    a.wall_v0[0,0,:] = np.array( [ 0., .0, 0. ] )
    a.wall_v0[0,1,:] = np.array( [ 0., -.0, 0. ] )
    a.wall_temp[0,0] = 0.33
    a.wall_temp[0,1] = 0.33
    a.init_r()
    a.init_v(.33)
    a.BC[0] = 1
    return a

def eight_PBC_wallx_tempgrad_8():
    """
    Returns a 8**3 system with density 10 and temperature T=.33.
    Periodic boundary conditions in y and z.
    Thermostatted walls in direction x with temperature of .13 and .53.
    """
    a = MPCD_system( (32,8,8) , 10, 1.)
    a.null_shift()
    print a
    a.tau = .5
    a.wall_v0[0,0,:] = np.array( [ 0., .0, 0. ] )
    a.wall_v0[0,1,:] = np.array( [ 0., -.0, 0. ] )
    a.wall_temp[0,0] = 0.37
    a.wall_temp[0,1] = 0.33
    a.init_r()
    a.init_v(.35)
    a.BC[0] = 1
    return a

def eight_PBC_wallx_sheary():
    """
    Returns a 8**3 system with density 10 and temperature T=.33.
    Periodic boundary conditions in y and z.
    Thermostatted walls in direction x applying a shear in the y direction.
    """
    a = MPCD_system( (8,8,8) , 10, 1.)
    a.null_shift()
    print a
    a.tau = .5
    a.wall_v0[0,0,:] = np.array( [ 0., .2, 0. ] )
    a.wall_v0[0,1,:] = np.array( [ 0., -.2, 0. ] )
    a.wall_temp[0,0] = 0.33
    a.wall_temp[0,1] = 0.33
    a.init_r()
    a.init_v(.33)
    a.BC[0] = 1
    return a

