
"""
MPCD.py - This file contains the base class to perform MPCD simulations: MPCD_system.
"""
## \namespace pyMPCD::MPCD
# \brief MPCD.py - This file contains the base class to perform MPCD simulations: MPCD_system.
#
# MPCD methods define the basic 

import numpy as np

from MPCD_f import mpcd_mod

## Defines all variables for a MPCD system
#
# A simulation box is either declared via MPCD_system or through one test_case or from a file.
# 
# \b Example
# \code
#import pyMPCD                               # import the pyMPCD package
#box = pyMPCD.MPCD_system( (8,8,8), 10, 1.0) # Create a PBC system of 8x8x8 cells, density 10 and cell size 1.0
# \endcode
class MPCD_system():
    """
    The MPCD_system class contains all variables relevant for MPCD simulations: box information, periodic boundaries, particles positions, ...
    """
    ## Initializes a MPCD_system
    #
    # \param N_cells a 3 elements list or tuple that contains the number of cells in each dimension.
    # \param density an integer number that define the average number of particles per cell.
    # \param a the linear cell size.
    def __init__( self, N_cells , density , a ):
        """
        Defines a MPCD_system with periodic boundary conditions.
        N_cells is the number of cells in the 3 dimensions.
        density is the reference density for initialization and for filling cells with virtual particles.
        a is the linear cell size.
        """
        ## Number of physical cells for the simulation.
        self.N_cells = np.array( N_cells , dtype=np.int32)
        # Check the input for N_cells
        if (len(self.N_cells) != 3): raise Exception
        ## The number of actual binning cells. Is higher than N_cells to allow for non-periodic systems.
        self.N_grid = self.N_cells + 1
        ## The average density (number of particles per cell).
        self.density = int(density)
        ## The total number of MPCD particles.
        self.so_N = np.prod( self.N_cells ) * self.density
        ## The linear cell size.
        self.a = float(a)
        ## The shift applied to the system.
        self.shift = np.zeros( (3,) , dtype=np.float64 )
        ## NumPy array for the position of the MPCD solvent.
        self.so_r = np.zeros( (self.so_N , 3) , dtype=np.float64 )
        ## A view to so_r in Fortran order.
        self.so_r_f = self.so_r.T
        ## NumPy array for the velocity of the MPCD solvent.
        self.so_v = np.zeros( (self.so_N , 3) , dtype=np.float64 )
        ## A view to so_v in Fortran order.
        self.so_v_f = self.so_v.T
        ## To delete
        self.r0 = np.zeros( (3,) , dtype=np.float64 )
        ## The size of the box.
        self.L = self.N_cells*self.a
        ## The MPCD time step, used for streaming.
        self.tau = float(0.)
        ## The number of particles in each cell.
        self.cells = np.zeros( self.N_grid, dtype=np.int32 )
        ## A view to cells in Fortran order.
        self.cells_f = self.cells.T
        self.par_list = np.zeros( (self.N_grid[0], self.N_grid[1], self.N_grid[2], 64), dtype=np.int32 )
        self.par_list_f = self.par_list.T
        self.v_com = np.zeros( (self.N_grid[0], self.N_grid[1], self.N_grid[2], 3), dtype=np.float64 )
        self.shift = np.zeros( (3,), dtype=np.float64)
        self.root = np.zeros( (3,), dtype=np.float64)
        self.BC = np.zeros( (3,) , dtype=np.int32 ) # type of wall. 0 = PBC , 1 = elastic collision with virtual particles.
        self.wall_v0 = np.zeros( (3, 2, 3) , dtype=np.float64 ) # mean velocity on each wall wall. indices = wall dir (x,y,z) , wall low/high , v
        self.wall_temp = np.zeros( (3, 2) , dtype=np.float64 ) # temperatures for the virtual particles on each wall side. indices = wall dir (x,y,z), wall low/high
        self.gravity = float(0.)

    def __str__(self):
        return str(type(self))+' size '+str(self.N_cells)+' , '+str(self.so_N)+' solvent particles'

    ## Initializes the particles according to a normal or flat velocity profile.
    # \param temp Initial temperature of the system.
    # \param boltz if True or unset, use a normal velocity profile of temperature T.
    #              Else, use a flat profile.
    def init_v(self, temp, boltz=True):
        """
        Initializes the particles according to a normal distribution of 
        temperature temp and resets the total velocity of the system.
        If boltz is set to False, a uniform distribution is used instead
        """
        if boltz:
            self.so_v[:,:] = np.random.randn( self.so_v.shape[0], self.so_v.shape[1] ) * np.sqrt(temp)
        else:
            self.so_v[:,:] = np.random.rand( self.so_v.shape[0], self.so_v.shape[1] )
            self.so_v[:,:] -= 0.5
            self.so_v[:,:] *= 2.*np.sqrt(6.*temp/2.)
        tot_v = np.sum( self.so_v , axis=0 ) / self.so_v.shape[0]
        tot_v = tot_v.reshape( ( 1 , tot_v.shape[0] ) )
        self.so_v -= tot_v

    ## Places particles in the simulation box at random according to a uniform distribution.
    def init_r(self):
        """
        Places particles in the simulation box at random according to a uniform distribution.
        """
        self.so_r[:,:] = np.random.rand( self.so_r.shape[0], self.so_r.shape[1] )
        self.so_r *= (self.a*self.N_cells).reshape( ( 1 , self.so_r.shape[1] ) )

    ## Advances the particles according to their velocities.
    def stream(self):
        """Advances the particles according to their velocities."""
        self.so_r[:] += self.so_v*self.tau

    ## Advances the particles according to their velocities and a constant acceleration in the z direction.
    # Also updates the velocities to take into account the acceleration.
    def accel(self):
        """
        Advances the particles according to their velocities and a constant acceleration in the z direction.
        Also updates the velocities to take into account the acceleration.
        """
        self.so_r[:] += self.so_v*self.tau + np.array( [0., 0., self.gravity] ).reshape( (1, 3) )*0.5*self.tau**2
        self.so_v[:] += np.array( [0., 0., self.gravity] ).reshape( (1, 3) )*self.tau


    def boundaries(self):
        """
        Corrects particles positions and velocities to take into account the boundary conditions.

        PBC keep the particles in the box by sending them to their periodic in-box location.
        Elastic walls reflect particles and reverse the velocities.
        """
        for i in range(3): # take each dim
            if (self.BC[i] == 0): # if PBC, simply appy x = mod( x , L ) to keep particles in the box
                self.so_r[:,i] = np.remainder( self.so_r[:,i] , self.L[i] )
            elif (self.BC[i] == 1): # if elastic wall, reflect "too high" particles around L and "too low" particles around 0
                j_out = ( self.so_r[:,i] > self.L[i] )
                self.so_r[j_out,i] = 2.*self.L[i] - self.so_r[j_out,i]
                self.so_v[j_out,i] = - self.so_v[j_out,i]
                j_out = ( self.so_r[:,i] < 0 )
                self.so_r[j_out,i] =  - self.so_r[j_out,i]
                self.so_v[j_out,i] =  - self.so_v[j_out,i]
            else:
                print "unknown boundary condition ", self.BC[i]
                raise Exception

    def check_in_box(self):
        """ 
        A test routine to check that all particles are actually in the box [0:L]
        """
        r_min = self.so_r.min(axis=0)
        t_min = ( r_min >= 0. ).min()
        r_max = self.so_r.max(axis=0)
        t_max = ( r_max < self.L ).min()
        if ( t_min and t_max ):
            return True
        else:
            return False
    def print_check_in_box(self):
        if (self.check_in_box()):
            print "All particles inside the box"
        else:
            print "Some particles outside the box"

    def null_shift(self):
        """
        Resets the shift to zero.
        """
        self.shift[:] = 0.
        self.root[:] = self.shift[:] - self.a

    def rand_shift(self):
        """
        Applies a random shift in [0:a[ to the system.
        """
        self.shift[:] = np.random.random( self.shift.shape[0] )*self.a
        self.root[:] = self.shift[:] - self.a

    def idx(self, i , cijk):
        """
        Returns in cijk the three cell indices for particle i.
        """
        np.floor( (self.so_r[i,:] - self.root) / self.a , cijk )
        for j in range(3):
            my_n = self.N_cells[j]
            if (self.BC[j] == 0):
                if ( cijk[j] >= my_n ): cijk[j] -= my_n

    def fill_box(self):
        """
        Bins the particles into cells.
        """
        cijk = np.zeros( (3,) , dtype=np.int32 )
        self.cells[:] = 0
        self.par_list[:] = 0
        for i in range(self.so_N):
            self.idx( i , cijk )
            my_n = self.cells[ cijk[0] , cijk[1] , cijk[2] ]
            self.par_list[ cijk[0] , cijk[1] , cijk[2] , my_n ] = i
            self.cells[ cijk[0] , cijk[1] , cijk[2] ] = my_n + 1
            
    def compute_v_com(self):
        """
        Computes the c.o.m. velocity for all cells.
        """
        self.v_com[:] = 0
        for ci in range(self.N_grid[0]):
            for cj in range(self.N_grid[1]):
                for ck in range(self.N_grid[2]):
                    v_local = np.zeros( (3, ) , dtype=np.float64)
                    n_local = self.cells[ci,cj,ck]
                    for i in range( n_local ):
                        part = self.par_list[ci,cj,ck,i]
                        v_local += self.so_v[part,:]
                    if (n_local > 0): self.v_com[ci,cj,ck,:] = v_local/n_local
        
    def MPCD_step_axis(self):
        """
        Performs a MPCD collision step where the axis is one of x, y or z, chosen at random.
        """
        v_therm = np.zeros((3,))
        nn = self.N_cells.copy()
        nn += self.BC
        is_wall = False
        v_wall = np.zeros( (3,) )
        for ci in range(nn[0]):
            for cj in range(nn[1]):
                for ck in range(nn[2]):
                    # Choose an axis to perform the rotation
                    rand_axis = np.random.randint(0,3)
                    axis1 = ( rand_axis + 1 ) % 3
                    axis2 = ( rand_axis + 2 ) % 3
                    if (np.random.randint(2)==0):
                        r_sign = 1
                    else:
                        r_sign = -1
                    # test if is a wall
                    is_wall = False
                    local_i = [ ci , cj , ck ]
                    for i in range(3):
                        if (self.BC[i]==1):
                            if ( (local_i[i]==0) or (local_i[i]==nn[i]-1) ):
                                is_wall=True
                                v_wall[:] = self.wall_v0[ i , min( local_i[i], 1 ) , : ]
                                v_temp = float(self.wall_temp[ i , min( local_i[i] , 1 ) ])
                    # number of particles in the cell
                    local_n = self.cells[ci,cj,ck]
                    # c.o.m. velocity in the cell
                    local_v = self.v_com[ci,cj,ck,:].copy()
                    # if cell is a wall, add virtual particles
                    if (is_wall):
                        if (local_n < self.density):
                            local_v = (
                                (np.random.randn(3) * np.sqrt(v_temp * (self.density - local_n) ) )
                                + v_wall*(self.density - local_n)
                                + local_v * local_n
                                ) / self.density
                            v_therm +=local_v
                    # perform cell-wise collisions
                    for i in range(local_n):
                        part = self.par_list[ci,cj,ck,i]
                        self.so_v[part,:] -= local_v
                        temp = self.so_v[part,axis2]
                        self.so_v[part,axis2] = r_sign*self.so_v[part,axis1]
                        self.so_v[part,axis1] = -r_sign*temp
                        self.so_v[part,:] += local_v
     
    ## Exchanges the positions and momenta of two solvent particles.
    # \param i index of the first particle to be exchanged.
    # \param j index of the second particle to be exchanged.
    def exchange_solvent(self,i,j):
        """
        Exchanges the positions and momenta of two solvent particles.
        """
        tmp_copy = self.so_r[i,:].copy()
        self.so_r[i,:] = self.so_r[j,:]
        self.so_r[j,:] = tmp_copy
        tmp_copy = self.so_v[i,:].copy()
        self.so_v[i,:] = self.so_v[j,:]
        self.so_v[j,:] = tmp_copy
        
    ## Sorts the solvent in the x,y,z cell order.
    def sort_solvent(self):
        """
        Sorts the solvent in the x,y,z cell order.
        """
        nn = self.N_cells.copy()
        nn += self.BC
        array_idx = 0
        for ci in range(nn[0]):
            for cj in range(nn[1]):
                for ck in range(nn[2]):
                    local_n = self.cells[ci,cj,ck]
                    for i in range(local_n):
                        self.exchange_solvent(self.par_list[ci,cj,ck,i],array_idx)
                        array_idx += 1

    def one_full_step(self):
        """
        Performs a full step of MPCD without gravitation, including the 
        streaming, taking into account the boundary conditions and the MPCD 
        collision step.
        """
        self.stream()
        self.boundaries()
        self.rand_shift()
        self.fill_box()
        self.compute_v_com()
        self.MPCD_step_axis()
    def one_full_accel(self):
        """
        Performs a full step of MPCD with gravitation, including the 
        streaming, taking into account the boundary conditions and the MPCD 
        collision step.
        """
        self.accel()
        self.boundaries()
        self.rand_shift()
        self.fill_box()
        self.compute_v_com()
        self.MPCD_step_axis()

    def one_full_step_f(self):
        """
        Performs a full step of MPCD without gravitation, including the 
        streaming, taking into account the boundary conditions and the MPCD 
        collision step.
        The streaming and binning steps are performed in Fortran.
        """
        #self.stream()
        mpcd_mod.stream(self.so_r_f, self.so_v_f, self.tau)
        self.boundaries()
        self.rand_shift()
        #self.fill_box()
        mpcd_mod.fill_box(self.so_r_f, self.cells_f, self.par_list_f, self.a, self.root)
        self.compute_v_com()
        self.MPCD_step_axis()

    


