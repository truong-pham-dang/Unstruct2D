! Author:  DANG Truong
! Purpose: to export grid in vtk format
! Date:    27/05/2017
subroutine visualize_grid
  use ModFiles
  use ModGeometry
  use ModInterfaces, only : DummyNodes, ErrorMessage
  implicit none

  integer :: i

  open(unit = 1,file = 'grid.vtk', status = 'replace')
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



  101   format('# vtk DataFile Version 2.0')
  102   format('VTK Format for unstructured grid')
  103   format('ASCII')
  104   format('DATASET UNSTRUCTURED_GRID')

  201   format('POINTS',I9,' float')
  202   format(3F25.12)
  203   format('CELLS',2I9)
  204   format (4I9)
  205   format('CELL_TYPES',I9)
  206   format(I9)

  close(1)

  return

end subroutine visualize_grid

! Author:  DANG Truong
! Purpose: to export grid in tecplot format
! Date:    04/06/2017
subroutine visualize_grid_tecplot
  use ModFiles
  use ModGeometry
  use ModInterfaces, only : DummyNodes, ErrorMessage
  implicit none

  integer :: i

  open(unit =2, file = 'mesh_tecplot.dat', status = 'replace', form = 'formatted')
  write(2,*) 'VARIABLES = "X", "Y"'
  write(2,*) 'ZONE NODES=',nndint,',ELEMENTS=',ntria, &
                  ',DATAPACKING = POINT, ZONETYPE = FETRIANGLE'

  do i = 1, nndint
    write(2,*) x(i),y(i)
  enddo

  do i = 1, ntria
   write(2,*) tria(1,i), tria(2,i), tria(3,i)
  enddo


  close(2)

end subroutine visualize_grid_tecplot
