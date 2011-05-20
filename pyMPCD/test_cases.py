from pyMPCD import MPCD_system
import numpy as np

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
    a = MPCD_system( (8,8,8) , 10, 1.)
    a.null_shift()
    print a
    a.tau = .5
    a.wall_v0[0,0,:] = np.array( [ 0., .0, 0. ] )
    a.wall_v0[0,1,:] = np.array( [ 0., -.0, 0. ] )
    a.wall_temp[0,0] = 0.33-2.
    a.wall_temp[0,1] = 0.33+2.
    a.init_r()
    a.init_v(.33)
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

