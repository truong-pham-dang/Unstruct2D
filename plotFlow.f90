!> @file plotFlow.f90
!!
!! Output of the flow field.
!
! *****************************************************************************
!
!  (c) J. Blazek, CFD Consulting & Analysis, www.cfd-ca.de
!  Created February 25, 2014
!  Last modification: June 4, 2014
!
! *****************************************************************************
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!
! *****************************************************************************

!> Writes out selected quantities for the whole flow field in Vis2D format.
!!
subroutine PlotFlow

  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none

! local variables
  character(chrlen) :: fname
  integer     :: errFlag, nquant, i, m
  real(rtype) :: rrho, u, v, e, press, temp, c, mach, ttot, ptot, machis, &
                 ptloss, pratio, ptotinf, gam1, ggm1, visc
  real(rtype) :: varout(mxqfield+2)

  character(20) :: citer
  integer     :: j_tecplot

! *****************************************************************************

  ptotinf = 0.D0
  if (kflow == "E") then
    gam1    = gamma - 1.D0
    ggm1    = gamma/gam1
    ptotinf = pinf*(1.D0+0.5D0*gam1*machinf*machinf)**ggm1
  endif

! open plot file

  write(fname,"(A,I5.5,A)") Trim(fnFlow),iter,".v2d"
  open(unit=ifFlow, file=fname, status="unknown", action="write", iostat=errFlag)
  if (errFlag /= 0) call ErrorMessage( "cannot open plot file (field)" )

! header

  nquant = 0
  do m=1,mxqfield
    if (lquant(m) == "Y") nquant = nquant + 1
  enddo

  write(ifFlow,1000) Trim(title),nquant+2

! names of variables

  do m=1,mxqfield
    if (lquant(m) == "Y") then
      write(ifFlow,"(A)") Trim(cquant(m))
    endif
  enddo

! number of data points and grid elements

  write(ifFlow,1010) nndint,ntria

! compute quantities & write'em out

  do i=1,nndint

    varout(1) = x(i)
    varout(2) = y(i)

    rrho  = 1.D0/cv(1,i)
    u     = cv(2,i)*rrho
    v     = cv(3,i)*rrho
    e     = cv(4,i)*rrho
    press = dv(1,i)
    temp  = dv(2,i)
    c     = dv(3,i)
    gam1  = dv(4,i) - 1.D0
    ggm1  = dv(4,i)/gam1

    if (kequs == "N") then
      visc = dv(6,i)
    else
      visc = 0.D0
    endif

    mach  = Sqrt(u*u+v*v)/c
    ttot  = (e+press*rrho)/dv(5,i)
    ptot  = press*(ttot/temp)**ggm1
    if (kflow == "E") then
      ptloss = 1.D0 - ptot/ptotinf
      pratio = ptotinf/press
    else
      ptloss = 1.D0 - ptot/ptinl
      pratio = ptinl/press
    endif
    machis = (pratio**(1.D0/ggm1)-1.D0)*2.D0/gam1
    machis = Max(machis, 0.D0)
    machis = Sqrt(machis)

! - store quantities in varout()
    nquant = 2
    do m=1,mxqfield
      if (lquant(m) == "Y") then
        nquant = nquant + 1

! ----- density
        if (m == 1) then
          varout(nquant) = cv(1,i)

! ----- u-velocity
        else if (m == 2) then
          varout(nquant) = u

! ----- v-velocity
        else if (m == 3) then
          varout(nquant) = v

! ----- static pressure
        else if (m == 4) then
          varout(nquant) = press

! ----- total pressure
        else if (m == 5) then
          varout(nquant) = ptot

! ----- static temperature
        else if (m == 6) then
          varout(nquant) = temp

! ----- total temperature
        else if (m == 7) then
          varout(nquant) = ttot

! ----- local Mach number
        else if (m == 8) then
          varout(nquant) = mach

! ----- isentropic Mach number
        else if (m == 9) then
          varout(nquant) = machis

! ----- total pressure loss
        else if (m == 10) then
          varout(nquant) = ptloss

! ----- laminar viscosity coefficient
        else if (m == 11) then
          varout(nquant) = visc
        endif
      endif
    enddo

    write(ifFlow,1020) (varout(m), m=1,nquant)

  enddo

! write node indexes

  do i=1,ntria
    write(ifFlow,1030) tria(1,i)-1,tria(2,i)-1,tria(3,i)-1
  enddo

  close(unit=ifFlow)

! Truong ajoute
!--------------------------------------------------------- 
  write(citer,'(I8.8)') iter 
  open(unit = 1,file = 'solution_'//trim(citer)//'.vtk', status = 'replace')
  write(1,101)
  write(1,102)
  write(1,103)
  write(1,104)
  
  write(1,201) nndint
  do i = 1, nndint
    write(1,202) x(i),y(i),0.0d0
  enddo

  write(1,203) ntria, 4*ntria

  do i = 1, ntria
	 write(1,204) 3, tria(1,i)-1,tria(2,i)-1,tria(3,i)-1
  enddo

  write(1,205) ntria
  do i = 1, ntria
    write(1,206) 5
  enddo

  write(1,302) nndint
! Density  
  write(1,301)
  write(1,300)
  do i = 1, nndint
    write(1,207) cv(1,i)
  enddo

! u-velocity
  write(1,303)
  write(1,300)
  do i = 1, nndint
    write(1,207) cv(2,i)/cv(1,i)
  enddo

! v-velocity
  write(1,304)
  write(1,300)
  do i = 1, nndint
    write(1,207) cv(3,i)/cv(1,i)
  enddo

! Velocity magnitude
  write(1,305)
  write(1,300)
  do i = 1, nndint
    write(1,207) sqrt((cv(2,i)/cv(1,i))**2+(cv(3,i)/cv(1,i))**2)
  enddo

! Pressure
  write(1,306)
  write(1,300)
  do i = 1, nndint
    write(1,207) dv(1,i)
  enddo
      
! Temperature
  write(1,307)
  write(1,300)
  do i = 1, nndint
    write(1,207) dv(2,i)
  enddo  

  close(1)

  open(unit = 2, file = 'solution_tecplot_'//trim(citer)//'.dat', status = 'replace')
  write(2,*) 'VARIABLES=X,Y,RHO,U,V,VELOCITY_MAGNITUDE,P,T'
  write(2,*) 'ZONE T="UNSTRUCTURED-COUNTOUR"'
  write(2,*) 'ZONETYPE=FEPOLYGON'
  write(2,*) 'NODES=',nndint
  write(2,*) 'ELEMENTS=',ntria
  write(2,*) 'FACES=',ntria*3
  write(2,*) 'NumConnectedBoundaryFaces=0'
  write(2,*) 'TotalNumBoundaryConnections=0'

  do i = 1, nndint
    write(2,401) x(i)
  enddo
  
  
  do i = 1, nndint
    write(2,401) y(i)
  enddo

! Density  
  do i = 1, nndint
    write(2,401) cv(1,i)
  enddo

! u-velocity
  do i = 1, nndint
    write(2,401) cv(2,i)/cv(1,i)
  enddo

! v-velocity
  do i = 1, nndint
    write(2,401) cv(3,i)/cv(1,i)
  enddo

! Velocity magnitude
  do i = 1, nndint
    write(2,401) sqrt((cv(2,i)/cv(1,i))**2+(cv(3,i)/cv(1,i))**2)
  enddo

! Pressure
  do i = 1, nndint
    write(2,401) dv(1,i)
  enddo
      
! Temperature
  do i = 1, nndint
    write(2,401) dv(2,i)
  enddo 
  
! Node indexes
  do i = 1, ntria
    write(2,*) tria(1,i), tria(2,i)
	write(2,*) tria(2,i), tria(3,i)
	write(2,*) tria(3,i), tria(1,i)
  enddo
  
  do i = 1, ntria
    write(2,204) i,i,i
  enddo  
  
  do i = 1, ntria
    write(2,204) 0,0,0
  enddo	 


  close(2) 

!----------------------------------------------------------

  

  101   format('# vtk DataFile Version 2.0')
  102   format('VTK Format for unstructed grid')
  103   format('ASCII')
  104   format('DATASET UNSTRUCTURED_GRID')

  201   format('POINTS',I9,' float')
  202   format(3F25.12)
  203   format('CELLS',2I9)
  204   format (4I9)
  205   format('CELL_TYPES',I9)
  206   format(I9)
  207   format(F25.18)

  300   format('LOOKUP_TABLE default')
  301   format('SCALARS RHO double')
  302   format('POINT_DATA',I9)
  303   format('SCALARS U double')
  304   format('SCALARS V double')
  305   format('SCALARS Velocity_magnitude double')
  306   format('SCALARS P double')
  307   format('SCALARS T double')

  401   format(F30.15)
    

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  format(1P,20E16.8)
1030  format(3I8)

end subroutine PlotFlow
