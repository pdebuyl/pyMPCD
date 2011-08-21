! Copyright (C) 2011 Pierre de Buyl

! This file is part of pyMPCD

! pyMPCD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! pyMPCD is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with pyMPCD.  If not, see <http://www.gnu.org/licenses/>.

!> This module provides drop-in Fortran replacement for some functions of the pyMPCD
!! package.
!!
!! The subroutines are given all variables explicitly, so that the code does not
!! depend on the setting of global variables in the Fortran module.

module mpcd_mod
  implicit none

contains
  
  !> Advances the particles according to their velocities.
  !!
  !! @param r fortran-ordered view of the positions of the MPCD particles.
  !! @param v fortran-ordered view of the velocities of the MPCD particles.
  !! @param tau timestep for the MPCD streaming.
  !! @param n number of particles of the r and v arrays. Becomes optional 
  !! in the f2py generated fortran object.
  subroutine stream(r, v, tau, n)
    implicit none
    double precision, intent(inout) :: r(3,n), v(3,n)
    double precision, intent(in) :: tau
    integer, intent(in) :: n

    integer :: i

    do i=1,n
       r(:,i) = r(:,i) + tau*v(:,i)
    end do

  end subroutine stream

  !> Fills the cells of a MPCD box with the indices of the MPCD particles
  !! that belong in them.
  !!
  !! @param r fortran-ordered view of the positions of the MPCD particles.
  !! @param cells number of particles in each MPCD cell.
  !! @param par_list the 0-based list of particles in each MPCD cell.
  !! @param a linear cell size.
  !! @param root the root of the 0th cell, including shifting.
  !! @param n number of particles of the r array. Becomes optional 
  !! in the f2py generated fortran object.
  !! @param Nx number of cells in the x-direction.
  !! @param Ny number of cells in the y-direction.
  !! @param Nz number of cells in the z-direction.
  subroutine fill_box(r, cells, par_list, a, root, n, Nx, Ny, Nz)
    double precision, intent(in) :: r(3,N)
    integer, intent(inout) :: cells(Nz,Ny,Nx)
    integer, intent(inout) :: par_list(64,Nz,Ny,Nx)
    double precision, intent(in) :: a, root(3)
    integer, intent(in) :: n
    integer, intent(in) :: Nx, Ny, Nz

    integer :: ci,cj,ck
    integer :: i

    cells = 0
    par_list = 0

    do i=1,n
       
       ci = floor( (r(1,i) - root(1)) / a ) + 1
       cj = floor( (r(2,i) - root(2)) / a ) + 1
       ck = floor( (r(3,i) - root(3)) / a ) + 1

       if ( ci > Nx .or. cj > Ny .or. ck > Nz ) then
          write(*,*) 'particle ', i, 'out of bounds'
       else
          cells(ck, cj, ci) = cells(ck, cj, ci) + 1
          if (cells(ck, cj, ci) > 64) then
             write(*,*) 'cell ',ci, cj,ck, ' full'
          else
             par_list(cells(ck,cj,ci), ck, cj, ci) = i-1
          end if
       end if

    end do


  end subroutine fill_box

  !>Exchanges the positions and velocities of two solvent particles.
  !! @param r array of positions.
  !! @param v array of velocities.
  !! @param species array of species.
  !! @param n number of particles in the arrays.
  !! @param i index of the first particle to be exchanged.
  !! @param j index of the second particle to be exchanged.
  subroutine exchange_solvent(r, v, species, n, i, j)
    implicit none
    double precision, intent(inout) :: r(3,N)
    double precision, intent(inout) :: v(3,N)
    integer, intent(inout) :: species(N)
    integer, intent(in) :: n
    integer, intent(in) :: i,j

    double precision :: x(3)
    integer :: s

    x = r(:,i)
    r(:,i) = r(:,j)
    r(:,j) = x

    x = v(:,i)
    v(:,i) = v(:,j)
    v(:,j) = x

    s = species(i)
    species(i) = species(j)
    species(j) = s

  end subroutine exchange_solvent

  !> Sorts the solvent in the x,y,z cell order.
  !! @param r fortran-ordered view of the positions of the MPCD particles.
  !! @param v fortran-ordered view of the velocities of the MPCD particles.
  !! @param species view of the species data of the MPCD particles.
  !! @param cells number of particles in each MPCD cell.
  !! @param par_list the 0-based list of particles in each MPCD cell.
  !! @param n number of particles of the r array. Becomes optional 
  !! in the f2py generated fortran object.
  !! @param Nx number of cells in the x-direction.
  !! @param Ny number of cells in the y-direction.
  !! @param Nz number of cells in the z-direction.
  subroutine sort_solvent(r, v, species, cells, par_list, n, Nx, Ny, Nz)
    implicit none
    double precision, intent(inout) :: r(3,N), v(3,N)
    integer, intent(inout) :: species(N)
    integer, intent(in) :: cells(Nz,Ny,Nx)
    integer, intent(in) :: par_list(64,Nz,Ny,Nx)
    integer, intent(in) :: n
    integer, intent(in) :: Nx, Ny, Nz

    integer :: ci,cj,ck,i,array_idx, j

    array_idx = 1
    do ci=1,Nx
       do cj=1,Ny
          do ck=1,Nz
             do i=1,cells(ck,cj,ci)
                j = par_list(i,ck,cj,ci)+1
                call exchange_solvent(r, v, species, n, j, array_idx)
                array_idx = array_idx+1
             end do
          end do
       end do
    end do

  end subroutine sort_solvent

end module mpcd_mod

