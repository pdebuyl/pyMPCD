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
MPCD.py - This file contains the base class to perform MPCD simulations: MPCD_system.
"""
## \namespace pyMPCD::MPCD
# \brief MPCD.py - This file contains the base class to perform MPCD simulations: MPCD_system.
#
# MPCD methods define the basic 

import numpy as np
from scipy.misc import factorial

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
    # \param N_species The number of solvent species to consider.
    def __init__( self, N_cells , density , a , N_species = 1):
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
        ## The number of solvent species
        self.N_species = N_species
        # Check that N_species is valid
        if (N_species < 1): raise Exception
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
        ## The local temperature in the cells.
        self.cells_temperature = np.zeros( (self.N_grid[0], self.N_grid[1], self.N_grid[2]), dtype=np.float64 )
        ## NumPy array for the species.
        self.so_species = np.zeros( (self.so_N, ), dtype=np.int32)
        ## NumPy array holding the mass of each species.
        self.so_mass = np.ones( (self.N_species, ), dtype=np.float64)
        ## The size of the box.
        self.L = self.N_cells*self.a
        ## The MPCD time step, used for streaming.
        self.tau = float(0.)
        ## The number of particles in each cell.
        self.cells = np.zeros( self.N_grid, dtype=np.int32 )
        ## A view to cells in Fortran order.
        self.cells_f = self.cells.T
        ## The list of particles, per cell.
        self.par_list = np.zeros( (self.N_grid[0], self.N_grid[1], self.N_grid[2], 64), dtype=np.int32 )
        ## A view to par_list in Fortran order.
        self.par_list_f = self.par_list.T
        ## The cell-wise center-of-mass velocity.
        self.v_com = np.zeros( (self.N_grid[0], self.N_grid[1], self.N_grid[2], 3), dtype=np.float64 )
        ## The origin of the grid.
        self.root = np.zeros( (3,), dtype=np.float64)
        ## Boundary conditions, in the three directions. 0 = PBC , 1 = elastic collision with virtual particles.
        self.BC = np.zeros( (3,) , dtype=np.int32 )
        ## Wall velocity, on each wall. indices = wall dir (x,y,z) , wall low/high , v.
        self.wall_v0 = np.zeros( (3, 2, 3) , dtype=np.float64 )
        ## Wall temperature for the virtual particles on each wall side. indices = wall dir (x,y,z), wall low/high.
        self.wall_temp = np.zeros( (3, 2) , dtype=np.float64 )
        ## Magnitude of the acceleration provided by gravity, if applicable.
        self.gravity = float(0.)
        ## Number of chemical reactions to consider.
        self.N_reactions = 0
        ## Kind of chemical reaction.
        self.reaction_kind = []
        ## Rates for the chemical reactions.
        self.rates = []
        ## Stoechiometry for reactants.
        self.reactants = []
        ## Stoechiometry for products.
        self.products = []

    def __str__(self):
        return str(type(self))+' size '+str(self.N_cells)+' , '+str(self.so_N)+' solvent particles'

    def add_reaction(self, rate, reactants, products):
        if ( not (len(reactants) == len(products) == self.N_species) ):
            print("Bad input for add_reaction")
            return
        if ( (min(reactants)<0) or (min(products)<0) ):
            print("Bad input for add_reaction")
            return
        N_reactants = np.sum(reactants)
        N_products = np.sum(products)
        if ( (N_reactants==1) and (N_products==1) ):
            kind = 0
        elif ( reactants[0]==1 and reactants[1]==2 and products[1]==3):
            kind = 1
        else:
            print("reaction type not supported in add_reaction")
            return

        self.N_reactions += 1
        self.reaction_kind.append(kind)
        self.rates.append(rate)
        self.reactants.append(reactants)
        self.products.append(products)

    def close_reaction(self):
        if (self.N_reactions > 0):
            self.rates = np.array(self.rates)
            self.reactants = np.array(self.reactants,ndmin=2)
            self.products = np.array(self.products,ndmin=2)

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
            self.so_v /= np.sqrt(self.so_mass[self.so_species]).reshape( (self.so_N, 1) )
        else:
            self.so_v[:,:] = np.random.rand( self.so_v.shape[0], self.so_v.shape[1] )
            self.so_v[:,:] -= 0.5
            self.so_v[:,:] *= 2.*np.sqrt(6.*temp/2.)
        tot_v = np.sum( self.so_v*(self.so_mass[self.so_species]).reshape( (self.so_N, 1) ) , axis=0 ) / self.so_mass[self.so_species].sum()
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

    ## Advances the particles according to their velocities by calling a Fortran routine.
    def stream_f(self):
        """"Advances the particles according to their velocities by calling a Fortran routine."""
        mpcd_mod.stream(self.so_r_f, self.so_v_f, self.tau)

    ## Advances the particles according to their velocities and a constant acceleration in the z direction.
    # Also updates the velocities to take into account the acceleration.
    def accel(self):
        """
        Advances the particles according to their velocities and a constant acceleration in the z direction.
        Also updates the velocities to take into account the acceleration.
        """
        self.so_r[:] += self.so_v*self.tau + np.array( [0., 0., self.gravity] ).reshape( (1, 3) )*0.5*self.tau**2
        self.so_v[:] += np.array( [0., 0., self.gravity] ).reshape( (1, 3) )*self.tau


    ## Corrects particles positions and velocities to take into account the boundary conditions.
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

    ## Bins the particles into the MPCD cells.
    def fill_box(self):
        """
        Bins the particles into the MPCD cells.
        """
        cijk = np.zeros( (3,) , dtype=np.int32 )
        self.cells[:] = 0
        self.par_list[:] = 0
        for i in range(self.so_N):
            self.idx( i , cijk )
            my_n = self.cells[ cijk[0] , cijk[1] , cijk[2] ]
            self.par_list[ cijk[0] , cijk[1] , cijk[2] , my_n ] = i
            self.cells[ cijk[0] , cijk[1] , cijk[2] ] = my_n + 1

    ## Bins the particles into the MPCD cells by calling a Fortran routine.
    def fill_box_f(self):
        """"Bins the particles into the MPCD cells by calling a Fortran routine."""
        mpcd_mod.fill_box(self.so_r_f, self.cells_f, self.par_list_f, self.a, self.root)

    ## Computes the center of mass velocity for all the cells in the system.
    # \param self A MPCD_system instance.
    def compute_v_com(self):
        """
        Computes the c.o.m. velocity for all cells.
        """
        self.v_com[:] = 0
        for ci in range(self.N_grid[0]):
            for cj in range(self.N_grid[1]):
                for ck in range(self.N_grid[2]):
                    mass_local = 0.
                    v_local = np.zeros( (3, ) , dtype=np.float64)
                    n_local = self.cells[ci,cj,ck]
                    for i in range( n_local ):
                        part = self.par_list[ci,cj,ck,i]
                        v_local += self.so_v[part,:]*self.so_mass[self.so_species[part]]
                        mass_local += self.so_mass[self.so_species[part]]
                    if (n_local > 0): self.v_com[ci,cj,ck,:] = v_local/mass_local


    ## Computes the temperature of the cells.
    # The temperature is computed as \f$ \frac{ \sum_{i=1}^{N^\xi} m_i (v_i-v_0) ^2 }{N^\xi-1} \f$
    # where \f$ v_0 \f$ is the center of mass velocity of the cell.
    # \param self A MPCD_system instance.
    def compute_cells_temperature(self):
        for ci in range(self.N_grid[0]):
            for cj in range(self.N_grid[1]):
                for ck in range(self.N_grid[2]):
                    v_local = self.v_com[ci,cj,ck,:]
                    T_local = 0.
                    n_local = self.cells[ci,cj,ck]
                    for i in range(n_local):
                        part = self.par_list[ci,cj,ck,i]
                        T_local += self.so_mass[self.so_species[part]]*((self.so_v[part,:]-v_local)**2).sum()
                    if (n_local>1):
                        self.cells_temperature[ci,cj,ck] = T_local/(3.*(n_local-1))
                    elif (n_local==1):
                        self.cells_temperature[ci,cj,ck] = 0.

    def MPCD_step_axis(self):
        """
        Performs a MPCD collision step where the axis is one of x, y or z, chosen at random.
        """
        v_therm = np.zeros((3,))
        nn = self.N_cells.copy()
        nn += self.BC
        is_wall = False
        v_wall = np.zeros( (3,) )
        local_N = np.zeros( (self.N_species,) , dtype=np.int32)
        amu = np.zeros( (self.N_reactions,), dtype=np.float64)
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
                    local_N *= 0
                    for i in range(local_n):
                        part = self.par_list[ci,cj,ck,i]
                        self.so_v[part,:] -= local_v
                        temp = self.so_v[part,axis2]
                        self.so_v[part,axis2] = r_sign*self.so_v[part,axis1]
                        self.so_v[part,axis1] = -r_sign*temp
                        self.so_v[part,:] += local_v
                        local_N[self.so_species[part]] += 1
                    # evaluate reaction probability
                    a0 = 0.
                    for i in range(self.N_reactions):
                        amu[i] = self.combi( i, local_N )*self.rates[i]
                        a0 += amu[i]
                    P_something = a0*self.tau
                    reac = -1
                    if (np.random.rand() < P_something):
                        amu /= a0
                        amu = np.cumsum(amu)
                        x = np.random.rand()
                        for i in range(self.N_reactions):
                            if (x < amu[i]):
                                reac = i
                                break
                    # apply reaction reac
                    if (reac>=0):
                        if (self.reaction_kind[reac]==0):
                            reac_s = np.where(self.reactants[i]==1)[0][0]
                            prod_s = np.where(self.products[i]==1)[0][0]
                            for i in range(local_n):
                                part = self.par_list[ci,cj,ck,i]
                                if (self.so_species[part]==reac_s):
                                    self.so_species[part] = prod_s
                                    break
                        elif (self.reaction_kind[reac]==1):
                            reac_s = np.where(self.reactants[i]==1)[0][0]
                            prod_s = np.where(self.products[i]==3)[0][0]
                            for i in range(local_n):
                                part = self.par_list[ci,cj,ck,i]
                                if (self.so_species[part]==reac_s):
                                    self.so_species[part] = prod_s
                                    break
                        else:
                            raise Exception

    def combi(self, reac, N):
        R = self.reactants[reac]
        P = self.products[reac]
        r = 1.0
        for i in range(self.N_species):
            if ( N[i] < R[i] ):
                return 0
            r *= factorial(N[i],exact=1)
            r /= factorial(N[i]-R[i],exact=1)
        return r
     
    ## Exchanges the positions, momenta and species of two solvent particles.
    # \param i index of the first particle to be exchanged.
    # \param j index of the second particle to be exchanged.
    def exchange_solvent(self,i,j):
        """
        Exchanges the positions, momenta and species of two solvent particles.
        """
        tmp_copy = self.so_r[i,:].copy()
        self.so_r[i,:] = self.so_r[j,:]
        self.so_r[j,:] = tmp_copy
        tmp_copy = self.so_v[i,:].copy()
        self.so_v[i,:] = self.so_v[j,:]
        self.so_v[j,:] = tmp_copy
        tmp_copy = self.so_species[i].copy()
        self.so_species[i] = self.so_species[j]
        self.so_species[j] = tmp_copy
        
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
        self.stream_f()
        self.boundaries()
        self.rand_shift()
        self.fill_box_f()
        self.compute_v_com()
        self.MPCD_step_axis()

    


