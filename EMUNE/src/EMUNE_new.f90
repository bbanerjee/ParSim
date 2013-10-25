!****************************************************************************
!*                                      *
!*  PROGRAM: EMUNE (ORIGINATED FROM 'EMU')                  *
!*      (2D DYNAMIC FRACTURE ANALYSIS CODE WITH BOND-BASED PERIDYNAMICS)*
!*                                      *
!*  PURPOSE: Given an Initial Grid, this program will generate a solution  *
!*      for this grid and based on a decision criteria, it will refine  *
!*      this grid and evaluate a new solution until the error < limit  *
!*                                      *
!*  PROGRAMMER: Dr. Youn Doh Ha                        *
!*                                      *
!****************************************************************************
!
!  EMUNE_new.f90 
!
!  FUNCTIONS:
!  EMUNE      - Entry point of console application.
!

program EMUNE

! use Main_Refine
  use Main_Solver
  use precision
  use Global_variables
  use Objects 
  use Input_subroutines

  real :: time1, time2

  character(len=25), parameter :: versn = 'EMUNE ver 1.0d 07/21/2009'
      
  call cpu_time(time1)
  open(6, file='emune.out',status='unknown') 
  write (6,*) versn  ! version information

  cont=0

! Read the main input file and the geometry of the body (nodes and elements)
  call read_input_file('main_input_file.txt')

! Perform the peridynamic analysis. It is the core of the program.
  call dynamic_solver        
  
  deallocate(nodes)
  deallocate(original_nodes)
  deallocate(elements)
  deallocate(quadtree_elements)
  deallocate(info_element)
  deallocate(original_quadtree_elements)

  call cpu_time(time2)
  print *, 'computational time(seconds):',time2-time1
  write(6,*) 'computational time(seconds):',time2-time1
  close(6)

end program EMUNE

