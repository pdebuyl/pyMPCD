module MPCD
  implicit none

contains

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

  subroutine fill_box(r, cells, par_list, a, root, n, Nx, Ny, Nz)
    double precision, intent(in) :: r(3,N)
    integer, intent(inout) :: cells(Nz,Ny,Nx)
    integer, intent(inout) :: par_list(64,Nz,Ny,Nx)
    double precision, intent(in) :: a, root(3)
    integer, intent(in) :: N
    integer, intent(in) :: Nx, Ny, Nz

    integer :: ci,cj,ck
    integer :: i

    cells = 0
    par_list = 0

    do i=1,N
       
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

end module MPCD

