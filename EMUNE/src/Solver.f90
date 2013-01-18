module Main_Solver    

  use Global_variables
  use dynamic_integration
  use precision
  use volume_partition
  use Horizon  
  use Objects
  use Input_subroutines

! Modification Summary (Mijia Yang)
 
! 09/12/06 - emu_main:  microcomputer version
 
! *** EMU ***

contains


  subroutine dynamic_solver

     integer(4) :: i,j,a,b
     character(len=25), parameter :: versn = 'EMU ver 1.0d 09/12/2006'
     character(len=100) :: mkdir_cmd, output_file_head;
     real(8) :: vel_top,vel_botton,time
     real(8) :: time3, time4
     logical :: flag

     ! Create a folder for the output data
     output_folder_name = './nodes_volume_output'
     mkdir_cmd = 'mkdir -p '//trim(output_folder_name)
     
     call execute_command_line(mkdir_cmd)

     ! Create a file name for the output data
     output_file_head = 'nodes_volume_output'
     output_file_name = trim(output_folder_name)//'/'//trim(output_file_head)
    
!    open(6, file='emu.out',status='unknown') 
    
!    open(69,file = 'velocity.txt', status='unknown')  
     open(666,file = 'nodes_volume_output_final.plt', status='unknown')
!    open(667,file = 'nnodes_in_iterations.txt', status='unknown')

!     write(666,*) 'TITLE="simulation results"'
!     write(666,*) 'VARIABLES="X","Y","DX","DY","VX","VY","DAM","W"'
!    write(666,*) 'VARIABLES="X","Y","DX","DY","VX","VY","DAM"'
    
!    open(11, file='history.out',status='unknown')
!    open(12, file='damage.txt', status='unknown')
!    open(15, file='damindex.txt', status='unknown')
!    open(10, file = 'emu_2d.txt')
    
    allocate(node_type(nnodes*mnode ))
    allocate(nofail(nnodes*mnode ))
!    allocate(omitt(nnodes*1 ))
    allocate(damage_index(nnodes*mnode ))

    allocate(fcoefm(nnodes*mnode ))
    allocate(fnorm(nnodes*mnode ))

!    allocate(broke(nnodes*1 ,nnodes*1 ))
    allocate(broke(maxfam ,nnodes*mnode))
  
    allocate(wt(nnodes*mnode ))
    allocate(f1(nnodes*mnode ))
    allocate(f2(nnodes*mnode ))

!    allocate(k1_x(nnodes*1))
!    allocate(k1_y(nnodes*1))
!    allocate(k2_x(nnodes*1))
!    allocate(k2_y(nnodes*1))
!    allocate(k3_x(nnodes*1))
!    allocate(k3_y(nnodes*1))
!    allocate(k4_x(nnodes*1))
!    allocate(k4_y(nnodes*1))
    
    allocate(edt_add(nnodes*mnode ))
    allocate(edt(nnodes*mnode ))
    allocate(vis(nnodes*mnode ))

!    allocate(family(nnodes*1 ))
    allocate(family(maxfam ))
    allocate(nodes_in_bin(nnodes*mnode ))

    allocate(radnod(nnodes*mnode ))
!    allocate(def_family(nnodes*1 ))
    allocate(def_family(maxfam ))
    allocate(nodes_in_bin_def(nnodes*mnode ))

      allocate(spsum(nnodes*mnode )) 
      allocate(ext1(nnodes*mnode ), ext2(nnodes*mnode )) 
      allocate(RadHorizon(nnodes*mnode )) 
      allocate(interval(nnodes*mnode, 2 )) 
      allocate(critical_strain(nnodes*mnode)) 

    allocate (xtempbotton(nnodes*mnode),ytempbotton(nnodes*mnode))
    allocate (xtemptop(nnodes*mnode), ytemptop(nnodes*mnode))
    allocate (ncorrespond_botton(nnodes*mnode), ncorrespond_top(nnodes*mnode))
    
      allocate(xt1(nnodes*mnode ), xt2(nnodes*mnode ), yt1(nnodes*mnode ), yt2(nnodes*mnode )) 
!    allocate(top_crack_speed(nt/snapshots_frequence +1))
!    allocate(botton_crack_speed(nt/snapshots_frequence +1))  

    do i=1,nnodes
      nodes(i)%disp(1)=zero
      nodes(i)%disp(2)=zero
      nodes(i)%new_disp(1)=zero
      nodes(i)%new_disp(2)=zero
      nodes(i)%veloc(1)=zero
      nodes(i)%veloc(2)=zero
      nodes(i)%new_veloc(1)=zero
      nodes(i)%new_veloc(2)=zero
      nodes(i)%accel(1) = zero
      nodes(i)%accel(2) = zero
      nodes(i)%strain_energy=zero
    enddo
    do i=1,nnodes*mnode
      ext1(i)=zero
      ext2(i)=zero
      wt(i)=zero
      f1(i)=zero
      f2(i)=zero
      vis(i)=zero
      spsum(i)=zero
      do j=1, maxfam
      broke(j,i) = 0
      enddo
    enddo
    

    node_type=1
    nofail=0
    damage_index=zero

      ecnode_temp = zero

     shortrange_dist_fac_nom = 1.35d0
     shortrange_dist_fac_init = 0.9d0
     shortrange_force_fac = 15.d0

   if (do_adap_refinement == 0) then
       call cal_volume       ! Calculate the volume of each node of the initial uniform grid
     call cal_interval       ! Calculate the horizon of each node the the initial uniform grid
   elseif (do_adap_refinement == 1) then
     print*,'...Calculating_Horizon...'
       call cpu_time(time3)
     call Calculate_Horizon() 
     call cpu_time(time4)
     print *, 'computational time(seconds):',time4-time3
     print*,'...Calculating_Volume...'
       call cpu_time(time3)
     call Calculate_Volume() 
     call cpu_time(time4)
     print *, 'computational time(seconds):',time4-time3
   endif

!   omitt=.FALSE.
   broke = 0
     vis = zero
!     call sort_ref

!   if (introduce_precrack == 1) call initial_broken_bond !>>03062009_YounDoh

!     call ExtForceDensity
!     call ExtForceDensity2
!   call InitialConditions

     call dynamics
     close(666)
!     close(6)
!     close(10)
!     close(11)
!   close(12)
!     close(15)


    deallocate(node_type, nofail)
    deallocate(damage_index, fcoefm, fnorm, broke)
    deallocate(wt, f1, f2)
    deallocate(edt_add, edt, vis)
    deallocate(family, nodes_in_bin, def_family, nodes_in_bin_def)

    deallocate(radnod, spsum, ext1, ext2)
      deallocate(RadHorizon, interval, critical_strain)

    deallocate(ncorrespond_botton, ncorrespond_top)
    deallocate(xtempbotton, ytempbotton, xtemptop, ytemptop)


 end subroutine dynamic_solver




End module Main_Solver    
