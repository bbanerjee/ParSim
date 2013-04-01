!********************************************************************
!********************************************************************
! ModuleFile : dynamic_integration
! Purpose    : This module includes subroutines to calculate the 
!              internal and external force density, which are used to 
!              calculate displacements
!********************************************************************
module dynamic_integration

  use Global_variables
  use Objects
  use Input_subroutines
! use Main_Refine
  use damage_computations
  use bond_family_computations
  use peridynamic_computations
  use precision
     
  use volume_partition
  use Horizon  
    
  implicit none

contains

  subroutine dynamics
  
    integer(4):: i
    ! integer(4):: j
    real :: time1, time2
    integer:: valid_node_count, output_file_count = 0
    character(len=100):: current_output_file_name

    ! omitt=.FALSE.
   
    ! do i=1, nnodes
    !   write(6,'(11F15.4)') iter*dt, nodes(i)%pos(1), nodes(i)%pos(2), nodes(i)%disp(1), nodes(i)%disp(2), nodes(i)%veloc(1), nodes(i)%veloc(2), nodes(i)%volume, nodes(i)%horizon_size, ext1(i), wt(i)
    ! enddo

    print*,"....Begin Solver...."

    do iter= 1, nt
      
      call cpu_time(time1)

      if (do_adap_refinement == 1) then ! If  is  =1 then the refinement and coarsening will be performed
        if(iter == 1) then
          call compute_critical_strain
          call sort_ref
          if (introduce_precrack == 1) call initial_broken_bond
!         call Refine
          broke = 0
          call compute_critical_strain
          call sort_ref
          if (introduce_precrack == 1) call initial_broken_bond
        endif
!       if((mod(iter,100) == 0).and.(iter>=40))then      
        if((mod(iter,10) == 0).and.(iter>=800))then      
!       if((mod(iter,5) == 0))then      
!         if(iter>1)then
!           call Coarse
!         endif
!         call Refine
          call compute_critical_strain
          call sort_ref
        endif
      else
        if (iter==1) then
          call compute_critical_strain
          call sort_ref
          if (introduce_precrack == 1) call initial_broken_bond
        endif
      endif

      call cpu_time(time2)
      print *, 'damage, computational time(seconds):',time2-time1

      ! Output nodal information every snapshots_freqence interation   
      if(iter==1)then
!       write(667,'(i6)') nnodes
!       do i=1,nnodes
!         write(666,'(3   f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%damage_index
!       enddo

        ! count valid nodes
        valid_node_count = 0
        do i = 1, nnodes
          if (.not.omitt(i)) then
            valid_node_count = valid_node_count + 1
          endif
        enddo

        write(666,*) 'TITLE="simulation results"'
        write(666,*) 'VARIABLES="X","Y","DX","DY","VX","VY","DAM","W"'
        write(666, '("ZONE I=",i5," F=POINT")') valid_node_count
        do i=1,nnodes
          if (.not.omitt(i)) then
!            write(666,'(7 f20.6)') nodes(i)%pos(1),nodes(i)%pos(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index
!            write(666,'(7 f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index
            write(666,'(8 f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index, wt(i)
          endif
        enddo
        write(*,*) 'snapshot at', iter
        write(6,*) 'snapshot at', iter

       
        ! Write the output to individual files
        write(current_output_file_name, fmt='(A,I5.5,A)') trim(output_file_name), output_file_count, '.tec'
        open(668, file=trim(current_output_file_name), status='unknown')
        write(668,*) 'TITLE="simulation results"'
        write(668,*) 'VARIABLES="X","Y","DX","DY","VX","VY","DAM","W"'
        write(668, fmt='(A,I5,A,F10.6,A)') 'ZONE I=',valid_node_count,' SOLUTIONTIME=', iter*dt,' F=POINT'
        do i=1,nnodes
          if (.not.omitt(i)) then
            write(668,'(8 f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index, wt(i)
          endif
        enddo
        close(668)
        output_file_count = output_file_count + 1

      endif  

      newiteration = .true.

      call cpu_time(time1)

      if ( iter == 1 ) then
!       ext1 = 0.d0
!       ext2 = 0.d0
!       call ExtForceDensity_Initial
!     elseif (iter == 2 ) then
        ext1 = zero
        ext2 = zero
        call ExtForceDensity
!       call ExtForceDensity_4_point_bending
      endif

      call cpu_time(time2)
      print *, 'ext, computational time(seconds):',time2-time1

      
!     call displace_boundary  ! Apply the displacement boundary conditions
   
!     Forward-Euler scheme
!     call peri_motion        ! Calculate the bond-forces in the grid and integrate the acceleration to obtain the new velocity of each node
!     call shortrange_motion  ! skip for crack_branching, don't skip for 4-point-bending
!     call intergrate_disp    ! integrate the new velocity to obtain the new displacement of each node

!     Velocity-Verlet scheme
      call cpu_time(time1)

!     call peri_motion_velocity_verlet
      call peri_motion_velocity_verlet2

      call cpu_time(time2)
      print *, 'peri_motion, computational time(seconds):',time2-time1

!     Runge-Kutta scheme
!     call peri_motion_RK     ! Calculate the bond-forces in the grid and integrate the acceleration to obtain the new velocity of each node
!     call intergrate_disp_RK ! integrate the new velocity to obtain the new displacement of each node
!     call intergrate_disp ! integrate the new velocity to obtain the new displacement of each node

      call boundary_cond      !>> 02262009_YounDoh
  
      ! Output nodal information every snapshots_freqence interation   
      if(mod(iter,snapshots_frequence) == 0)then
!       write(667,'(i6)') nnodes
!       do i=1,nnodes
!         write(666,'(3   f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%damage_index
!       enddo
        
        ! count valid nodes
        valid_node_count = 0
        do i = 1, nnodes
          if (.not.omitt(i)) then
            valid_node_count = valid_node_count + 1
          endif
        enddo

        write(666,*) 'TITLE="simulation results"'
        write(666,*) 'VARIABLES="X","Y","DX","DY","VX","VY","DAM","W"'
        write(666, '("ZONE I=",i5," F=POINT")') valid_node_count
        do i=1,nnodes
          if (.not.omitt(i)) then
!            write(666,'(7 f20.6)') nodes(i)%pos(1),nodes(i)%pos(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index
!            write(666,'(7 f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index
            write(666,'(8 f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index, wt(i)
          endif
        enddo
        write(*,*) 'snapshot at', iter
        write(6,*) 'snapshot at', iter

        ! Write the output to individual files
        write(current_output_file_name, fmt='(A,I5.5,A)') trim(output_file_name), output_file_count, '.tec'
        open(668, file=trim(current_output_file_name), status='unknown')
        write(668,*) 'TITLE="simulation results"'
        write(668,*) 'VARIABLES="X","Y","DX","DY","VX","VY","DAM","W"'
        write(668, fmt='(A,I5,A,F10.6,A)') 'ZONE I=',valid_node_count,' SOLUTIONTIME=', iter*dt,' F=POINT'
        do i=1,nnodes
          if (.not.omitt(i)) then
            write(668,'(8 f20.6)') nodes(i)%pos(1)+nodes(i)%disp(1),nodes(i)%pos(2)+nodes(i)%disp(2),nodes(i)%disp(1),nodes(i)%disp(2),nodes(i)%veloc(1),nodes(i)%veloc(2), nodes(i)%damage_index, wt(i)
          endif
        enddo
        close(668)
        output_file_count = output_file_count + 1
      endif  

        
    enddo

  end subroutine dynamics


!___02262009_____________________________________________    
  subroutine boundary_cond
    ! apply displacement or velocity boundary conditions.
    integer(4):: i, mbd
    real(8):: xt1, xt2, ut1, ut2, vt1, vt2, xt1new, xt2new, ut1new, ut2new, vt1new, vt2new

    mbd = 0  ! mbd == 1: apply boundary condition, mpd == 0: no B.C

    if (mbd>0) then
      do i=1,nnodes
        xt1    = nodes(i)%pos(1)
        xt2    = nodes(i)%pos(2)
        xt1new = nodes(i)%pos(1) + nodes(i)%disp(1)
        xt2new = nodes(i)%pos(2) + nodes(i)%disp(2)
        ut1    = nodes(i)%disp(1)
        ut2    = nodes(i)%disp(2)
        ut1new = nodes(i)%new_disp(1)
        ut2new = nodes(i)%new_disp(2)
        vt1    = nodes(i)%veloc(1)
        vt2    = nodes(i)%veloc(2)
        vt1new = nodes(i)%new_veloc(1)
        vt2new = nodes(i)%new_veloc(2)
        ! prescribed displacement (4-point-bending)
        if((abs(xt2)<=0.0042).and.(abs(xt1-0.397)<=0.0042))then
          nodes(i)%new_veloc(2) = zero
          nodes(i)%new_disp(2) = zero
        elseif((abs(xt2)<=0.0042).and.(abs(xt1-0.916)<=0.0042))then
          nodes(i)%new_veloc(2) = zero
          nodes(i)%new_disp(2) = zero
        endif
        ! prescribed velocity (4-point-bending)
!       if((abs(xt2-0.306).le.1.e-6).and.(abs(xt1-0.0).le.1.e-6))then
!         nodes(i)%new_veloc(2) = -bc_velocity*0.13
!!        if (iter.le.120000) then
!!          nodes(i)%new_veloc(2) = -bc_velocity*0.13*real(iter)/real(120000)
!!        else
!!          nodes(i)%new_veloc(2) = -bc_velocity*0.13
!!        endif
!         nodes(i)%new_disp(2) = nodes(i)%disp(2)+ nodes(i)%new_veloc(2)*dt
!       elseif((abs(xt2-0.306).le.1.e-6).and.(abs(xt1-0.519)<=0.005))then
!!        nodes(i)%new_veloc(2) = -bc_velocity*1.13
!         nodes(i)%new_veloc(2) = -bc_velocity
!!        if (iter.le.120000) then
!!          nodes(i)%new_veloc(2) = -bc_velocity*1.13*real(iter)/real(120000)
!!        else
!!          nodes(i)%new_veloc(2) = -bc_velocity*1.13
!!        endif
!         nodes(i)%new_disp(2) = nodes(i)%disp(2)+ nodes(i)%new_veloc(2)*dt
!       endif
        ! prescribed velocity (crack_branching)
!       if(abs(xt2).le.1.e-6) then ! Nodes at botton boundary 
!         nodes(i)%new_veloc(2) = -bc_velocity
!         nodes(i)%new_disp(2) = nodes(i)%disp(2)+ nodes(i)%new_veloc(2)*dt
!       elseif(abs(xt2-0.04).le.1.e-6) then ! Nodes at top boundary 
!         nodes(i)%new_veloc(2) = bc_velocity
!         nodes(i)%new_disp(2) = nodes(i)%disp(2)+ nodes(i)%new_veloc(2)*dt
!       endif
!       nodes(53)%new_veloc(1) = 0.d0
!       nodes(53)%new_disp(1) = 0.d0
!       nodes(1166)%new_veloc(1) = 0.d0
!       nodes(1166)%new_disp(1) = 0.d0

        ! prescribed displacement (simple stretching plate)
!       if((abs(xt2)<=1.e-4))then
!         nodes(i)%new_veloc(2) = 0.d0
!         nodes(i)%new_disp(2) = 0.d0
!       endif
!       if((abs(xt1)<=1.e-4))then
!         nodes(i)%new_veloc(1) = 0.d0
!         nodes(i)%new_disp(1) = 0.d0
!       endif
      end do
    end if

    ! Replace the values of displacement, velocity and strain energy with new ones
    do i=1,nnodes
      nodes(i)%disp(1)=nodes(i)%new_disp(1)
      nodes(i)%disp(2)=nodes(i)%new_disp(2)
      nodes(i)%veloc(1)=nodes(i)%new_veloc(1)
      nodes(i)%veloc(2)=nodes(i)%new_veloc(2)
      nodes(i)%strain_energy=wt(i)
    enddo
    
  end subroutine boundary_cond

!********************************************************************

  subroutine ExtForceDensity_Initial

    integer(4) :: i, j, icount1, icount2, num, num2
    ! real(8) :: xtempbotton(nnodes)
    real(8) :: ytempbotton(nnodes), xtemptop(nnodes), ytemptop(nnodes)
    integer(4) :: ncorrespond_botton(nnodes), ncorrespond_top(nnodes)
    real(8) :: ctran
    icount1=0
    icount2=0
    ctran=0.0
    do i=1,nnodes
!     if(abs(nodes(i)%pos(1)).le.1.0e-6) then ! Nodes at botton boundary
!       icount1=icount1+1
!       xtempbotton(icount1)=nodes(i)%pos(1)
!       ytempbotton(icount1)=nodes(i)%pos(2)
!       ncorrespond_botton(icount1)=i
!!      nofail(i) = 1
!     endif
!     if(abs(nodes(i)%pos(1)-0.8).le.1.0e-6) then ! Nodes at top boundary
      if(abs(nodes(i)%pos(1)-0.4).le.1.0d-6) then ! Nodes at top boundary
        icount2=icount2+1
        xtemptop(icount2)=nodes(i)%pos(1)
        ytemptop(icount2)=nodes(i)%pos(2)
        ncorrespond_top(icount2)=i
!       nofail(i) = 1
      endif
      nofail(i) = 1
    end do
!   do i=1, icount1
!     do j=i+1, icount1
!       if(xtempbotton(j).gt.xtempbotton(i)) then
!         ctran=xtempbotton(i)
!         xtempbotton(i)=xtempbotton(j)
!         xtempbotton(j)=ctran
!       endif
!     enddo
!   enddo
!   do i=1, icount2
!     do j=i+1, icount2
!       if(xtemptop(j).gt.xtemptop(i)) then
!         ctran=xtemptop(i)
!         xtemptop(i)=xtemptop(j)
!         xtemptop(j)=ctran
!       endif
!     enddo
!   enddo
    do i=1, icount1
      do j=i+1, icount1
        if(ytempbotton(j).gt.ytempbotton(i)) then
          ctran=ytempbotton(i)
          ytempbotton(i)=ytempbotton(j)
          ytempbotton(j)=ctran
        endif
      enddo
    enddo
    do i=1, icount2
      do j=i+1, icount2
        if(ytemptop(j).gt.ytemptop(i)) then
          ctran=ytemptop(i)
          ytemptop(i)=ytemptop(j)
          ytemptop(j)=ctran
        endif
      enddo
    enddo

    do i=1,icount1  ! vanished the thickness effect
      num=ncorrespond_botton(i)
      if(i.ge.2.and.i.le.icount1-1) then
        ext1(num)=-force_mag*abs(ytempbotton(i+1)-ytempbotton(i-1))/(nodes(num)%volume*2.0d0)

      else
        if(i == 1) ext1(num)=-force_mag*abs(ytempbotton(i+1)-ytempbotton(i))/(nodes(num)%volume*2.0d0)
        if(i == icount1) ext1(num)=-force_mag*abs(ytempbotton(i)-ytempbotton(i-1))/(nodes(num)%volume*2.0d0)
      endif
    enddo
     
    num =ncorrespond_top(1)
    num2=0
    do j=1,num-num2
      ctran  = force_mag*abs(ytemptop(2)-ytemptop(1))/(nodes(num)%volume*2.0d0)
      ext1(j+num2)= ctran * dble(j-1)/dble(num-num2-1)
    enddo
    do i=1,icount2
      num=ncorrespond_top(i)
      if(i.ge.2.and.i.le.icount2-1) then
        num2=ncorrespond_top(i-1)
        do j=1,num-num2
          ctran  = force_mag*abs(ytemptop(i+1)-ytemptop(i-1))/(nodes(num)%volume*2.0d0)
          ext1(j+num2)= ctran * dble(j-1)/dble(num-num2-1)
        enddo
!       ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i-1))/(nodes(num)%volume*2.0d0)
!     else
!       if(i == 1) ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i))/(nodes(num)%volume*2.0d0)
!       if(i == icount2) ext1(num)=force_mag*abs(ytemptop(i)-ytemptop(i-1))/(nodes(num)%volume*2.0d0)
      endif
    enddo
    num =ncorrespond_top(icount2)
    num2=ncorrespond_top(icount2-1)
    do j=1,num-num2
      ctran  = force_mag*abs(ytemptop(icount2)-ytemptop(icount2-1))/(nodes(num)%volume*2.0d0)
      ext1(j+num2)= ctran * dble(j-1)/dble(num-num2-1)
    enddo
    
  end subroutine ExtForceDensity_Initial

!********************************************************************
! subroutine ExtForceDenstiy
! Purpose : set the external force density array
! Define variables
!      ForceA -- define the location of load
!      b      -- external force density array
!********************************************************************
  subroutine ExtForceDensity

    integer(4) :: i, j, icount1, icount2, num
!   real(8) :: xtempbotton(nnodes),ytempbotton(nnodes), xtemptop(nnodes), ytemptop(nnodes)
!   real(8), dimension(:), allocatable :: xtempbotton(:), ytempbotton(:), xtemptop(:), ytemptop(:)
!   integer(4) :: ncorrespond_botton(nnodes), ncorrespond_top(nnodes)
!   integer(4), dimension(:), allocatable :: ncorrespond_botton(:), ncorrespond_top(:)
    real(8) :: ctran
  
!   allocate (xtempbotton(nnodes),ytempbotton(nnodes), xtemptop(nnodes), ytemptop(nnodes))
!   allocate (ncorrespond_botton(nnodes), ncorrespond_top(nnodes))
  
    icount1=0
    icount2=0
    ctran=0.0
    do i=1,nnodes
!     if(abs(nodes(i)%pos(1)).le.1.0e-6) then ! Nodes at botton boundary
!     if(abs(nodes(i)%pos(2)).le.1.0d-6) then
      if(abs(nodes(i)%pos(2)+0.02d0).le.1.0d-6) then
        icount1=icount1+1
        xtempbotton(icount1)=nodes(i)%pos(1)
        ytempbotton(icount1)=nodes(i)%pos(2)
        ncorrespond_botton(icount1)=i
        nofail(i) = 1
      endif
!     if(abs(nodes(i)%pos(1)-0.8).le.1.0e-6) then ! Nodes at top boundary
!     if(abs(nodes(i)%pos(1)-0.4).le.1.0e-6) then ! Nodes at top boundary
!     if(abs(nodes(i)%pos(2)-0.04d0).le.1.0d-6) then
      if(abs(nodes(i)%pos(2)-0.02d0).le.1.0d-6) then
        icount2=icount2+1
        xtemptop(icount2)=nodes(i)%pos(1)
        ytemptop(icount2)=nodes(i)%pos(2)
        ncorrespond_top(icount2)=i
        nofail(i) = 1
      endif
!     nofail(i) = 1
    end do
    do i=1, icount1
      do j=i+1, icount1
        if(xtempbotton(j).gt.xtempbotton(i)) then
          ctran=xtempbotton(i)
          xtempbotton(i)=xtempbotton(j)
          xtempbotton(j)=ctran
        endif
      enddo
    enddo
    do i=1, icount2
      do j=i+1, icount2
        if(xtemptop(j).gt.xtemptop(i)) then
          ctran=xtemptop(i)
          xtemptop(i)=xtemptop(j)
          xtemptop(j)=ctran
        endif
      enddo
    enddo
!   do i=1, icount1
!     do j=i+1, icount1
!       if(ytempbotton(j).gt.ytempbotton(i)) then
!         ctran=ytempbotton(i)
!         ytempbotton(i)=ytempbotton(j)
!         ytempbotton(j)=ctran
!       endif
!     enddo
!   enddo
!   do i=1, icount2
!     do j=i+1, icount2
!       if(ytemptop(j).gt.ytemptop(i)) then
!         ctran=ytemptop(i)
!         ytemptop(i)=ytemptop(j)
!         ytemptop(j)=ctran
!       endif
!     enddo
!   enddo

    do i=1,icount1  ! vanished the thickness effect
      num=ncorrespond_botton(i)
      if(i.ge.2.and.i.le.icount1-1) then
!       ext1(num)=-force_mag*abs(ytempbotton(i+1)-ytempbotton(i-1))/(nodes(num)%volume*2.0d0)
        ext2(num)=-force_mag*dabs(xtempbotton(i+1)-xtempbotton(i-1))/(nodes(num)%volume*2.0d0)
!       ext1(num)=-force_mag
      else
!       if(i == 1) ext1(num)=-force_mag*abs(ytempbotton(i+1)-ytempbotton(i))/(nodes(num)%volume*2.0d0)
!       if(i == icount1) ext1(num)=-force_mag*abs(ytempbotton(i)-ytempbotton(i-1))/(nodes(num)%volume*2.0d0)
        if(i == 1) ext2(num)=-force_mag*dabs(xtempbotton(i+1)-xtempbotton(i))/(nodes(num)%volume*2.d0)
        if(i == icount1) ext2(num)=-force_mag*dabs(xtempbotton(i)-xtempbotton(i-1))/(nodes(num)%volume*2.d0)
!       if(i == 1) ext1(num)=-force_mag
!       if(i == icount1) ext1(num)=-force_mag
      endif
    enddo
     
    do i=1,icount2
      num=ncorrespond_top(i)
      if(i.ge.2.and.i.le.icount2-1) then
!       ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i-1))/(nodes(num)%volume*2.0d0)
        ext2(num)=force_mag*dabs(xtemptop(i+1)-xtemptop(i-1))/(nodes(num)%volume*2.0d0)
!       ext1(num)=force_mag
      else
!       if(i == 1) ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i))/(nodes(num)%volume*2.0d0)
!       if(i == icount2) ext1(num)=force_mag*abs(ytemptop(i)-ytemptop(i-1))/(nodes(num)%volume*2.0d0)
!       if(i == 1) ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i))/(nodes(num)%volume*4.0d0)
!       if(i == icount2) ext1(num)=force_mag*abs(ytemptop(i)-ytemptop(i-1))/(nodes(num)%volume*4.0d0)
        if(i == 1) ext2(num)=force_mag*dabs(xtemptop(i+1)-xtemptop(i))/(nodes(num)%volume*2.0d0)
        if(i == icount2) ext2(num)=force_mag*dabs(xtemptop(i)-xtemptop(i-1))/(nodes(num)%volume*2.0d0)
!       if(i == 1) ext1(num)=force_mag
!       if(i == icount2) ext1(num)=force_mag
      endif
    enddo
  
!   deallocate(ncorrespond_botton, ncorrespond_top)
!   deallocate(xtempbotton, ytempbotton, xtemptop, ytemptop)
    
  end subroutine ExtForceDensity

  subroutine ExtForceDensity2

    integer(4) :: i, j, icount1, icount2, num
!   real(8) :: xtempbotton(nnodes),ytempbotton(nnodes), xtemptop(nnodes), ytemptop(nnodes)
!   integer(4) :: ncorrespond_botton(nnodes), ncorrespond_top(nnodes)
    real(8) :: ctran
    icount1=0
    icount2=0
    ctran=0.0
    do i=1,nnodes
      if(abs(nodes(i)%pos(1)).le.1.0e-6) then ! Nodes at left boundary
        icount1=icount1+1
        xtempbotton(icount1)=nodes(i)%pos(1)
        ytempbotton(icount1)=nodes(i)%pos(2)
        ncorrespond_botton(icount1)=i
        nofail(i) = 1
      endif
      if(abs(nodes(i)%pos(1)-0.1).le.1.0e-6) then ! Nodes at right boundary
        icount2=icount2+1
        xtemptop(icount2)=nodes(i)%pos(1)
        ytemptop(icount2)=nodes(i)%pos(2)
        ncorrespond_top(icount2)=i
        nofail(i) = 1
      endif
    end do
    do i=1, icount1
      do j=i+1, icount1
        if(ytempbotton(j).gt.ytempbotton(i)) then
          ctran=ytempbotton(i)
          ytempbotton(i)=ytempbotton(j)
          ytempbotton(j)=ctran
        endif
      enddo
    enddo
    do i=1, icount2
      do j=i+1, icount2
        if(ytemptop(j).gt.ytemptop(i)) then
          ctran=ytemptop(i)
          ytemptop(i)=ytemptop(j)
          ytemptop(j)=ctran
        endif
      enddo
    enddo

    do i=1,icount1  ! vanished the thickness effect
      num=ncorrespond_botton(i)
      if(i.ge.2.and.i.le.icount1-1) then
        ext1(num)=-force_mag*abs(ytempbotton(i+1)-ytempbotton(i-1))/2.0d0
      else
        if(i == 1) ext1(num)=-force_mag*abs(ytempbotton(i+1)-ytempbotton(i))/2.0d0
        if(i == icount1) ext1(num)=-force_mag*abs(ytempbotton(i)-ytempbotton(i-1))/2.0d0
      endif
    enddo
     
    do i=1,icount2
      num=ncorrespond_top(i)
      if(i.ge.2.and.i.le.icount2-1) then
        ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i-1))/2.0d0
      else
        if(i == 1) ext1(num)=force_mag*abs(ytemptop(i+1)-ytemptop(i))/2.0d0
        if(i == icount2) ext1(num)=force_mag*abs(ytemptop(i)-ytemptop(i-1))/2.0d0
      endif
    enddo
    
  end subroutine ExtForceDensity2

  subroutine ExtForceDensity_4_point_bending

    integer(4) :: i, j, icount1, icount2, num, ntran
!   real(8) :: xtempbotton(nnodes),ytempbotton(nnodes), xtemptop(nnodes), ytemptop(nnodes)
!   integer(4) :: ncorrespond_botton(nnodes), ncorrespond_top(nnodes)
    real(8) :: ctran, force_mag2, force_time
    icount1=0
    icount2=0
    ctran=zero
    ntran=0
    do i=1,nnodes
!     if((abs(nodes(i)%pos(2)-0.306).le.1.e-6).and.(abs(nodes(i)%pos(1)).le.0.0042)) then
      if((abs(nodes(i)%pos(2)-0.306d0).le.1.d-6).and.(abs(nodes(i)%pos(1)).le.0.0082d0)) then
!     if((abs(nodes(i)%pos(2)-0.306).le.1.e-6).and.(abs(nodes(i)%pos(1)).le.0.0122)) then
        icount1=icount1+1
        xtempbotton(icount1)=nodes(i)%pos(1)
        ytempbotton(icount1)=nodes(i)%pos(2)
        ncorrespond_botton(icount1)=i
!       nofail(i) = 1
      endif
!     if((abs(nodes(i)%pos(2)-0.306).le.1.e-6).and.(abs(nodes(i)%pos(1)-0.519).le.0.0042))then
      if((abs(nodes(i)%pos(2)-0.306d0).le.1.d-6).and.(abs(nodes(i)%pos(1)-0.519d0).le.0.0082d0))then
!     if((abs(nodes(i)%pos(2)-0.306).le.1.e-6).and.(abs(nodes(i)%pos(1)-0.519).le.0.0122))then
        icount2=icount2+1
        xtemptop(icount2)=nodes(i)%pos(1)
        ytemptop(icount2)=nodes(i)%pos(2)
        ncorrespond_top(icount2)=i
!       nofail(i) = 1
      endif
!     nofail(i) = 1
    end do
    do i=1, icount1
      do j=i+1, icount1
        if(xtempbotton(j).gt.xtempbotton(i)) then
          ctran=xtempbotton(i)
          xtempbotton(i)=xtempbotton(j)
          xtempbotton(j)=ctran
          ntran=ncorrespond_botton(i)
          ncorrespond_botton(i)=ncorrespond_botton(j)
          ncorrespond_botton(j)=ntran
        endif
      enddo
    enddo
    do i=1, icount2
      do j=i+1, icount2
        if(xtemptop(j).gt.xtemptop(i)) then
          ctran=xtemptop(i)
          xtemptop(i)=xtemptop(j)
          xtemptop(j)=ctran
          ntran=ncorrespond_top(i)
          ncorrespond_top(i)=ncorrespond_top(j)
          ncorrespond_top(j)=ntran
        endif
      enddo
    enddo

!   vanished the thickness effect
    do i=1,nnodes

      if(abs(nodes(i)%pos(1)).le.0.000001d0)then
        nofail(i) = 1
      endif
      if(abs(nodes(i)%pos(1)-0.916d0).le.0.000001d0)then
        nofail(i) = 1
      endif
!     if(abs(nodes(i)%pos(2)-0.306).le.0.005)then
!       nofail(i) = 1
!     endif
      if((abs(nodes(i)%pos(2)).le.0.000001d0))then
        nofail(i) = 1
      endif

!     if((abs(nodes(i)%pos(2)-0.306).le.0.000001).and.(abs(nodes(i)%pos(1)).le.0.000001))then
!!      nofail(i) = 1
!       ext2(i)=-0.13*force_mag*abs(nodes(i+1)%pos(1)-nodes(i)%pos(1))/(nodes(i)%volume*2.0d0)
!     elseif((abs(nodes(i)%pos(2)-0.306).le.0.000001).and.(abs(nodes(i)%pos(1)-0.519)<=0.005))then
!!      nofail(i) = 1
!       ext2(i)=-force_mag*abs(nodes(i+1)%pos(1)-nodes(i-1)%pos(1))/(nodes(i)%volume*2.0d0)
!     endif
    enddo

    force_time = dble(iter)*dt
!   if (force_time .le. real(90)*dt) then
!     force_mag2 = real(iter)/90.0*force_mag
!   else
      force_mag2 = force_mag
!   endif

    do i=1,icount1  ! vanished the thickness effect
      num=ncorrespond_botton(i)
      if(i.ge.2.and.i.le.icount1-1) then
        ext2(num)=-0.13*force_mag2*dabs(xtempbotton(i+1)-xtempbotton(i-1))/(nodes(num)%volume*2.0d0)
      else
        if(i == 1) ext2(num)=-0.13*force_mag2*dabs(xtempbotton(i+1)-xtempbotton(i))/(nodes(num)%volume*2.d0)
        if(i == icount1) ext2(num)=-0.13*force_mag2*dabs(xtempbotton(i)-xtempbotton(i-1))/(nodes(num)%volume*2.d0)
      endif
    enddo
    do i=1,icount2
      num=ncorrespond_top(i)
      if(i.ge.2.and.i.le.icount2-1) then
        ext2(num)=-force_mag2*dabs(xtemptop(i+1)-xtemptop(i-1))/(nodes(num)%volume*2.0d0)
      else
        if(i == 1) ext2(num)=-force_mag2*dabs(xtemptop(i+1)-xtemptop(i))/(nodes(num)%volume*2.0d0)
        if(i == icount2) ext2(num)=-force_mag2*dabs(xtemptop(i)-xtemptop(i-1))/(nodes(num)%volume*2.0d0)
      endif
!     print *, ext2(num)
    enddo
    
  end subroutine ExtForceDensity_4_point_bending

!********************************************************************

  !-- intergrate_disp
  subroutine intergrate_disp

    !real(8) visco >>03072009_YounDoh
    integer(4):: m


!   do m=1,nnodes
!     nodes(m)%new_veloc(1) = nodes(m)%new_veloc(1)*0.8d0
!     nodes(m)%new_veloc(2) = nodes(m)%new_veloc(2)*0.8d0
!   end do
!   visco=0.2  ! apply the damage viscosity (visco)
!   visco=0.0
    do m=1,nnodes
      nodes(m)%new_veloc(1) = nodes(m)%new_veloc(1)/(One+visco*vis(m))
      nodes(m)%new_veloc(2) = nodes(m)%new_veloc(2)/(One+visco*vis(m))
    end do
    ! integrate for new displacements
    do m=1,nnodes
      if(omitt(m)) then
        nodes(m)%new_veloc(1) = 0.
        nodes(m)%new_veloc(2) = 0.
        nodes(m)%new_disp(1) = nodes(m)%disp(1)
        nodes(m)%new_disp(2) = nodes(m)%disp(2)
      else
        nodes(m)%new_disp(1) = nodes(m)%disp(1)+nodes(m)%veloc(1)*dt
        nodes(m)%new_disp(2) = nodes(m)%disp(2)+nodes(m)%veloc(2)*dt
      end if
    end do
!>> move to boundary_cond
    ! Update the values of displacement, velocity and strain energy
!   do m=1,nnodes
!     nodes(m)%disp(1)=nodes(m)%new_disp(1)
!     nodes(m)%disp(2)=nodes(m)%new_disp(2)
!     nodes(m)%veloc(1)=nodes(m)%new_veloc(1)
!     nodes(m)%veloc(2)=nodes(m)%new_veloc(2)
!     nodes(m)%strain_energy=wt(m)
!   enddo
!<< 02262009_YounDoh
        
  end subroutine intergrate_disp

!********************************************************************
!********************************************************************
! subroutine InitialConditions
! Purpose : set initial condition of the problem on the grid
! Define variables
!  
!********************************************************************
  ! Initial Condition
  
  subroutine InitialConditions

    implicit none
    
    integer :: i, icount1, icount2
    ! integer :: j, num
!   real(8) :: xtempbotton(nnodes),ytempbotton(nnodes), xtemptop(nnodes), ytemptop(nnodes)
!   integer :: ncorrespond_botton(nnodes), ncorrespond_top(nnodes)
    real(8) :: ctran
    icount1=0
    icount2=0
    ctran=zero
    do i=1,nnodes
        
      nodes(i)%new_disp(1) = zero
      nodes(i)%disp(1) = zero   
      nodes(i)%veloc(1) = zero
      nodes(i)%new_veloc(1) = zero

      nodes(i)%new_disp(2) = zero
      nodes(i)%disp(2) = zero   

!     4-point-bending
!     if((abs(nodes(i)%pos(2)-0.306).le. 0.000001).and.(abs(nodes(i)%pos(1)-0.0).le. 0.000001))then
!       nodes(i)%veloc(2) = -bc_velocity*0.13
!     elseif((abs(nodes(i)%pos(2)-0.306).le. 0.000001).and.(abs(nodes(i)%pos(1)-0.519)<=0.005))then
!       nodes(i)%veloc(2) = -bc_velocity
!     else
!       nodes(i)%veloc(2) = 0.0d0
!     endif
!     nodes(i)%new_veloc(2) = 0.0d0

!     crack_branching   
!     nodes(i)%veloc(2) = bc_velocity*(nodes(i)%pos(2)-0.02)/0.02
!     nodes(i)%new_veloc(2) = 0.0d0
      nodes(i)%veloc(2) = zero
      nodes(i)%new_veloc(2) = zero

    end do

  end subroutine InitialConditions


!********************************************************************
 end module dynamic_integration
!********************************************************************
