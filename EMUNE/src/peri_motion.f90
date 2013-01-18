!****************************************************************************
!****************************************************************************
!*                                      *
!*  ModuleFile : peridynamic_computations                  *
!*  Purpose    : This module includes subroutines to calculate the      *
!*        internal and external force density, which are used to    *
!*        calculate displacements                    *
!*                                      *
!*  PROGRAMMER : Dr. Youn Doh Ha                      *
!*                                      *
!****************************************************************************
module peridynamic_computations

  use Global_variables
  use Objects
  use Input_subroutines
  use precision
  use bond_family_computations

  implicit integer(4) (i-n)
  implicit real(8) (a-h,o-z)

contains

subroutine peri_motion_velocity_verlet

  ! Integrates the node accelerations due to peridynamic ("structured") interaction.
  ! Update the node velocities
  ! Also update the bond offset strains for microplastic materials
  ! break bonds based on critical strian for breakage 

  n_broken_bounds = 0
  if( dynamic_or_static == 1) then  ! dynamic
    dampx =  zero
    dampy =  zero
  else                ! quasi-static (dynamic relaxation)
    dampx =  0.2d0
    dampy =  0.2d0
  endif

  twt = zero
  ! Loop over all the nodes.
  do mi=1,nnodes
    if(omitt(mi)) goto 100
    if(nodes(mi)%iflag .eqv. .TRUE.) goto 100

    u1i = nodes(mi)%disp(1)
    u2i = nodes(mi)%disp(2)
    x1i = nodes(mi)%pos(1)
    x2i = nodes(mi)%pos(2)
    v1i = nodes(mi)%veloc(1)
    v2i = nodes(mi)%veloc(2)
    u1io= nodes(mi)%old_disp(1)
    u2io= nodes(mi)%old_disp(2)
    nofli = nofail(mi)   !>>030509_YounDoh

    f1(mi) = zero
    f2(mi) = zero
    wt(mi) = zero
    delta=nodes(mi)%horizon_size

    ! Get the family of node mi (all the nodes within its horizon, delta).
    call get_family(mi)

    ! brokij is set to 1 if the bond between mi and mj is already broken.
    nbroke=0

    ! Loop over nodes in the family of node mi.
    do js=1,mfam
      mj = family(js)

      if(nodes(mj)%iflag.eqv. .TRUE.) continue  ! hanging node, no calculation

      nodtj = node_type(mj)
      
!      brokij = broke(mi,mj)
      brokij = broke(js,mi)
!      if (introduce_precrack==1.and.(.not.brokij))then
!        call broke_reg(mi,mj)
!      endif
        
      if((.not.omitt(mi)).and.(.not.omitt(mj)).and.(.not.brokij)) then
        u1j = nodes(mj)%disp(1)
        u2j = nodes(mj)%disp(2)
        x1j = nodes(mj)%pos(1)
        x2j = nodes(mj)%pos(2)
        v1j = nodes(mj)%veloc(1)
        v2j = nodes(mj)%veloc(2)
        u1jo= nodes(mj)%old_disp(1)
        u2jo= nodes(mj)%old_disp(2)
        vr1 = v1j-v1i
        vr2 = v2j-v2i
        noflj = nofail(mj)

        fcof=One 
        ! Find the peridynamic interparticle force.

        call peri_force(u1j-u1i,u2j-u2i,x1j-x1i,x2j-x2i,brokij,expco,df1,df2,dw,spcoef,str)

        ! Critical stretch as a function of damage.
        dmgij = dmax1(damage_index(mi),damage_index(mj))
!        coef1 = 0.2d0
!        coef2 = 0.2d0
!        coef3 = 1.4d0
        coef1 = dc1
        coef2 = dc2
        coef3 = dc3
        if(coef2>zero .and. coef3>One) then
          if(dmgij<=coef1) then
            fac = One
          elseif(dmgij<0.9999d0) then
            fac = One+coef2*(dmgij-coef1)/(One-dmgij)
            fac = dmin1(fac, coef3)
          else
            fac = coef3
          endif
        else
          fac = One
        endif

                
!        ecr2 = critical_stretch*fac
        ecr2 = critical_strain(mi)*fac
        ! Break bond if critical stretch exceeded.
!                   if(str>critical_stretch.and..not.brokij) then
        if(str>ecr2.and..not.brokij) then
          if(nofli==0.and.noflj==0) then
!            broke(mi,mj) = 1
            broke(js,mi) = 1
            vis(mi) = vis(mi)+One
            n_broken_bounds = n_broken_bounds+1
          endif
        end if

        ! Sum up the force on node mi.
        dvolj = nodes(mj)%volume

        ! reduce volume of node mj if it is not fully within the horizon.
        xi = dsqrt( (x1j-x1i)**2 + (x2j-x2i)**2 )
        radij = 0.5d0*dmax1(interval(mj,1),interval(mj,2))

        if(radij > zero) then
          if(xi <= delta-radij) then
            fac = One
          elseif(xi <= delta+radij) then
            fac = (delta+radij-xi)/(2.d0*radij)
          else
            fac = zero
          endif
        else
          fac = zero
        endif

        df1 = df1*dvolj*fac
        df2 = df2*dvolj*fac
        f1(mi) = f1(mi)+df1
        f2(mi) = f2(mi)+df2
        wt(mi) = wt(mi)+0.5d0*dw*dvolj*fac
        spsum(mi) = spsum(mi)+dvolj*spcoef/denst
        ! Sum up the force across the planes x1 = x1force, ...

      end if

!      if(broke(mi,mj).eqv. .TRUE.1) then
      if(broke(js,mi) == 1) then
        nbroke=nbroke+1
      endif
    end do

    if (mfam.gt.0) then
      f1b=dble(nbroke)
      f2b=dble(mfam)
      nodes(mi)%damage_index = f1b/f2b ! IF index is 0 then the node has no broken bond
      damage_index(mi)=f1b/f2b
    else
      nodes(mi)%damage_index = zero
      damage_index(mi) = zero
    endif
      
100    if(omitt(mi)) then
      nodes(mi)%new_veloc(1) = zero
      nodes(mi)%new_veloc(2) = zero
      nodes(mi)%new_disp(1) = nodes(mi)%disp(1)
      nodes(mi)%new_disp(2) = nodes(mi)%disp(2)
    else
!      fac = denst/dt
!      write(*,*) fac,nodes(mi)%veloc(1)
!      nodes(mi)%new_veloc(1) = (ext1(mi)+f1+fac*nodes(mi)%veloc(1))/(fac+dampx)
!      nodes(mi)%new_veloc(2) = (ext2(mi)+f2+fac*nodes(mi)%veloc(2))/(fac+dampy)
      if( dynamic_or_static == 1) then ! dynamic
        if (iter == 1) then ! Forward-Euler scheme for the first step
          fac = dt/denst
          nodes(mi)%new_disp(1) = nodes(mi)%disp(1)+nodes(mi)%veloc(1)*dt
          nodes(mi)%new_disp(2) = nodes(mi)%disp(2)+nodes(mi)%veloc(2)*dt
          nodes(mi)%new_veloc(1) = nodes(mi)%veloc(1)
          nodes(mi)%new_veloc(2) = nodes(mi)%veloc(2)
        else
          fac = dt*dt/denst
!            fac = dt*dt/denst*0.8d0
          nodes(mi)%new_disp(1) = 2.d0*nodes(mi)%disp(1)-nodes(mi)%old_disp(1)+(ext1(mi)+f1(mi))*fac
          nodes(mi)%new_disp(2) = 2.d0*nodes(mi)%disp(2)-nodes(mi)%old_disp(2)+(ext2(mi)+f2(mi))*fac
          nodes(mi)%new_veloc(1) = (nodes(mi)%new_disp(1)-nodes(mi)%old_disp(1))/(2.d0*dt)
          nodes(mi)%new_veloc(2) = (nodes(mi)%new_disp(2)-nodes(mi)%old_disp(2))/(2.d0*dt)
        endif
      endif
    end if

    nodes(mi)%old_disp(1) = nodes(mi)%disp(1)
    nodes(mi)%old_disp(2) = nodes(mi)%disp(2)

!    twt = twt + wt(mi)*nodes(mi)%volume
  end do
!  print *, twt

  ! find the viscosity coefficient.
  do mi=1,nnodes
    if(.not.omitt(mi)) vis(mi) = (One-visdk)*vis(mi)
  enddo
  
  return

end subroutine peri_motion_velocity_verlet
!_____________________________________________________________________________________
subroutine peri_motion_velocity_verlet2

  ! Integrates the node accelerations due to peridynamic ("structured") interaction.
  ! Update the node velocities
  ! Also update the bond offset strains for microplastic materials
  ! break bonds based on critical strian for breakage 
  ! [one-step Velocity-Verlet formulation]
  ! 1. v(n+1/2) = v(n) + dt/2m * f(q(n))
  ! 2. q(n+1) = q(n) + dt * v(n+1/2)
  ! 3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))

!    write (*,'(f25.12)') dt
!    write (*,'(f25.12)') young
!    write (*,'(f25.12)') denst
!    write (*,'(f25.12)') fracture_energy
!    write (*,'(f25.12)') horizon_factor
!    write (*,'(f25.12)') initial_crack_x0
!    write (*,'(f25.12)') initial_crack_y0
!    write (*,'(f25.12)') initial_crack_x1
!    write (*,'(f25.12)') initial_crack_y1
!    write (*,'(f25.12)') force_mag
!    write (*,'(f25.12)') dc3

  n_broken_bounds = 0
  if( dynamic_or_static == 1) then  ! dynamic
    dampx =  zero
    dampy =  zero
  else                ! quasi-static (dynamic relaxation)
    dampx =  0.2d0
    dampy =  0.2d0
  endif

    ! force update at previous configuration (n)
    if ((iter == 1).or.(modified_mesh .eqv. .TRUE.)) then
        do mi=1,nnodes
           ! Get the family of node mi (all the nodes within its horizon, delta).
          delta=nodes(mi)%horizon_size
           call get_family(mi)

          u1i = nodes(mi)%disp(1)
          u2i = nodes(mi)%disp(2)
          x1i = nodes(mi)%pos(1)
          x2i = nodes(mi)%pos(2)

          f1(mi) = zero
          f2(mi) = zero

           ! Loop over nodes in the family of node mi.
           do js=1,mfam
             mj = family(js)
             brokij = broke(js,mi)
            
             if((.not.omitt(mi)).and.(.not.omitt(mj)).and.(.not.brokij)) then
               u1j = nodes(mj)%disp(1)
              u2j = nodes(mj)%disp(2)
              x1j = nodes(mj)%pos(1)
              x2j = nodes(mj)%pos(2)

              ! Find the peridynamic interparticle force.
              call peri_force(u1j-u1i,u2j-u2i,x1j-x1i,x2j-x2i,brokij,expco,df1,df2,dw,spcoef,str)

               ! Sum up the force on node mi.
              dvolj = nodes(mj)%volume

              ! reduce volume of node mj if it is not fully within the horizon.
              xi = dsqrt( (x1j-x1i)**2 + (x2j-x2i)**2 )
              radij = 0.5d0*dmax1(interval(mj,1),interval(mj,2))

               if(radij > zero) then
                if(xi <= delta-radij) then
                  fac = One
                elseif(xi <= delta+radij) then
                  fac = (delta+radij-xi)/(2.d0*radij)
!                           fac = zero  ! without volume correction
                else
                   fac = zero
                endif
              else
                fac = zero
              endif

               df1 = df1*dvolj*fac
              df2 = df2*dvolj*fac
!                   force at the current configuration (n+1)
              f1(mi) = f1(mi)+df1
              f2(mi) = f2(mi)+df2

            end if
           end do
        end do
    endif

!   1st step & 2nd step (intermediate velocity update and position update)
    do mi=1,nnodes
         if(omitt(mi).or.(nodes(mi)%iflag.eqv. .TRUE.)) then
        nodes(mi)%new_veloc(1) = zero
        nodes(mi)%new_veloc(2) = zero
        nodes(mi)%new_disp(1) = nodes(mi)%disp(1)
          nodes(mi)%new_disp(2) = nodes(mi)%disp(2)
         else
          if( dynamic_or_static == 1) then ! dynamic
!               1. v(n+1/2) = v(n) + dt/2m * f(q(n))
            fac = dt/denst*0.5d0
             nodes(mi)%new_veloc(1) = nodes(mi)%veloc(1) + (ext1(mi)+f1(mi))*fac
            nodes(mi)%new_veloc(2) = nodes(mi)%veloc(2) + (ext2(mi)+f2(mi))*fac
                 nodes(mi)%veloc(1) = nodes(mi)%new_veloc(1)
                 nodes(mi)%veloc(2) = nodes(mi)%new_veloc(2)
!               2. q(n+1) = q(n) + dt * v(n+1/2)
            nodes(mi)%new_disp(1)  = nodes(mi)%disp(1)  + nodes(mi)%veloc(1)*dt
            nodes(mi)%new_disp(2)  = nodes(mi)%disp(2)  + nodes(mi)%veloc(2)*dt
           endif
        end if
    enddo


!   3rd step (force update and velocity update)
!   3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
!   twt = zero

    do mi=1,nnodes
      ! Get the family of node mi (all the nodes within its horizon, delta).
        delta=nodes(mi)%horizon_size
      call get_family(mi)

      u1i = nodes(mi)%new_disp(1)
      u2i = nodes(mi)%new_disp(2)
      x1i = nodes(mi)%pos(1)
      x2i = nodes(mi)%pos(2)
      nofli = nofail(mi)

      f1(mi) = zero
      f2(mi) = zero
      wt(mi) = zero
      spsum(mi) = zero

      ! brokij is set to 1 if the bond between mi and mj is already broken.
      nbroke=0

      ! Loop over nodes in the family of node mi.
      do js=1,mfam
        mj = family(js)
        brokij = broke(js,mi)

        if((.not.omitt(mi)).and.(.not.omitt(mj)).and.(.not.brokij)) then
          u1j = nodes(mj)%new_disp(1)
          u2j = nodes(mj)%new_disp(2)
          x1j = nodes(mj)%pos(1)
          x2j = nodes(mj)%pos(2)
          noflj = nofail(mj)

          ! Find the peridynamic interparticle force.
          call peri_force(u1j-u1i,u2j-u2i,x1j-x1i,x2j-x2i,brokij,expco,df1,df2,dw,spcoef,str)

          ! Critical stretch as a function of damage.
          dmgij = dmax1(damage_index(mi),damage_index(mj))
          coef1 = dc1
          coef2 = dc2
          coef3 = dc3
          if(coef2>zero .and. coef3>One) then
            if(dmgij<=coef1) then
              fac = One
            elseif(dmgij<0.9999d0) then
              fac = One+coef2*(dmgij-coef1)/(One-dmgij)
              fac = dmin1(fac, coef3)
            else
              fac = coef3
            endif
          else
            fac = One
          endif

          ! Break bond if critical stretch exceeded.
          ecr2 = critical_strain(mi)*fac
          if((str>ecr2).and.(.not.brokij)) then
            if(nofli==0.and.noflj==0) then
              broke(js,mi) = 1
              vis(mi) = vis(mi)+One
              n_broken_bounds = n_broken_bounds+1
            endif
          end if

          ! Sum up the force on node mi.
          dvolj = nodes(mj)%volume

          ! reduce volume of node mj if it is not fully within the horizon.
          xi = dsqrt( (x1j-x1i)**2 + (x2j-x2i)**2 )
          radij = 0.5d0*dmax1(interval(mj,1),interval(mj,2))

          if(radij > zero) then
            if(xi <= delta-radij) then
              fac = One
            elseif(xi <= delta+radij) then
              fac = (delta+radij-xi)/(2.d0*radij)
!                        fac = zero  ! without volume correction
            else
              fac = zero
            endif
          else
            fac = zero
          endif

          df1 = df1*dvolj*fac
          df2 = df2*dvolj*fac
!               force at the current configuration (n+1)
          f1(mi) = f1(mi)+df1
          f2(mi) = f2(mi)+df2
          wt(mi) = wt(mi)+0.5d0*dw*dvolj*fac
          spsum(mi) = spsum(mi)+dvolj*spcoef/denst
        end if

        if(broke(js,mi) == 1) then
          nbroke=nbroke+1
        endif
      end do

      if (mfam.gt.0) then
        f1b=dble(nbroke)
        f2b=dble(mfam)
        nodes(mi)%damage_index = f1b/f2b ! IF index is 0 then the node has no broken bond
        damage_index(mi)=f1b/f2b
      else
        nodes(mi)%damage_index = zero
        damage_index(mi) = zero
      endif
          
      if(omitt(mi).or.(nodes(mi)%iflag.eqv. .TRUE.)) then
           nodes(mi)%new_veloc(1) = zero
           nodes(mi)%new_veloc(2) = zero
           nodes(mi)%new_disp(1) = nodes(mi)%disp(1)
          nodes(mi)%new_disp(2) = nodes(mi)%disp(2)
      else
        if( dynamic_or_static == 1) then ! dynamic
          fac = dt/denst*0.5d0
          nodes(mi)%new_veloc(1) = nodes(mi)%veloc(1) + (ext1(mi)+f1(mi))*fac
          nodes(mi)%new_veloc(2) = nodes(mi)%veloc(2) + (ext2(mi)+f2(mi))*fac
        endif
      end if

!      twt = twt + wt(mi)*nodes(mi)%volume

    end do

!     print *, twt

  ! find the viscosity coefficient.
  do mi=1,nnodes
    if(.not.omitt(mi)) vis(mi) = (One-visdk)*vis(mi)
  enddo
  
  return

end subroutine peri_motion_velocity_verlet2

subroutine peri_force(eta1,eta2,xi1,xi2,brokp,expco,df1,df2,dw,spcoef,str)

  logical brokp 

  ! returns force density per unit volume due to peridynamic interaction between nodes

  ! r is the initial distance between the nodes
  r = dsqrt(xi1**2+xi2**2)
  p1 = xi1+eta1
  p2 = xi2+eta2

  ! p is the new distance between the nodes
  if(dsqrt(p1**2+p2**2).gt.1.0d16) then
    stop 1111   ! diverge
  endif
  p = dsqrt(p1**2+p2**2)

  call MicromodF(r,dmicroF, delta,young)

  ! elasticity model.
  ! find bond strain.
  ! find the force
  if(p>zero.and.r<delta) then
    u = p-r ! u is the displacement between the two nodes
    if(r>zero.and.r<delta) then
      str = (p-r)/r
    else
      str = zero
    end if
 
    df1 = dmicroF*str*p1/p
    df2 = dmicroF*str*p2/p
    dw = 0.5d0*dmicroF*u**2
    spcoef = dmicroF
!    linearized model
!    df1 = microF*(eta1*xi1*xi1/(r**2)+eta2*xi1*xi2/(r**2))
!    df2 = microF*(eta1*xi1*xi2/(r**2)+eta2*xi2*xi2/(r**2))
!    dw = 0.5d0*microF*u**2
!    spcoef = microF
  else
    df1 = zero
    df2 = zero
    dw = zero
    spcoef = zero
  end if

  return
end subroutine peri_force

subroutine MicromodF(esi,dmicroF, RadH,young)

! Purpose : define the micromodulus function C
! Define variables:
!  kesi   -- the relative position in reference configuration
!  microF -- the micromodulus function
!  real, parameter :: thickness=1.0
    
  if (sw_micro == 0) then
!  2D constant micromodulus   
!    dmicroF =  6.0d0*young/(pi*thickness*(RadH**3)*(1.d0/3.d0)*(4.d0/3.d0))  ==> thickness effect will be vanished in volume integration
!    dmicroF =  6.0d0*young/(pi*(RadH**3)*(1.d0/3.d0)*(4.d0/3.d0))  ==> simplify
!    dmicroF =  6.0d0*young/(pi*(RadH**3)*(4.d0/9.d0))  ==> simplify
    dmicroF =  13.5d0*young/(pi*(RadH**3))
  elseif (sw_micro == 1) then
!  2D conical micromodulus
!    dmicroF = 24.0d0*young*(1.0d0-esi/RadH)/(pi*thickness*(RadH**3)*(1.d0/3.d0)*(4.d0/3.d0))  ==> thickness effect will be vanished in volume integration
!    dmicroF = 24.0d0*young*(1.0d0-esi/RadH)/(pi*(RadH**3)*(1.d0/3.d0)*(4.d0/3.d0))  ==> simplify
!    dmicroF = 24.0d0*young*(1.0d0-esi/RadH)/(pi*(RadH**3)*(4.d0/9.d0))  ==> simplify
    dmicroF = 54.0d0*young*(1.0d0-esi/RadH)/(pi*(RadH**3))
  endif

end subroutine MicromodF


subroutine shortrange_motion
 
      ! Integrates the equation of motion for short-range forces ("structureless interaction") between nodes.
 
      ! Structureless interaction between target nodes:
      ! For each node on this processor, find the deformed family.
      ! Sum up the structureless force from each node in the deformed family.


      ! Sort the nodes on this processor according to deformed position.
      call sort_def()

      do mi=1,nnodes
         if(nodes(mi)%iflag.eqv. .TRUE.) continue

            x1i = nodes(mi)%pos(1)
            x2i = nodes(mi)%pos(2)
            y1i = nodes(mi)%pos(1)+nodes(mi)%disp(1)
            y2i = nodes(mi)%pos(2)+nodes(mi)%disp(2)
            v1i = nodes(mi)%veloc(1)
            v2i = nodes(mi)%veloc(2)          
            dvoli = nodes(mi)%volume
            nodti = node_type(mi)
            xmassi = denst*dvoli
!<<
!            radnodi = nodes(mi)%horizon_size/horizon_factor
            radnodi = radnod(mi)
!>>02272009_YounDoh

          ! Find the deformed family for node m. This is all the nodes that
          ! could possibly have a structureless interaction with node m.

          rad_search = 3.d0*radnod_max
!          rad_search = 3.0*radnodi
!          call get_family(mi)
          call get_def_family(mi)


          ! Loop over nodes in the deformed family.
!<<
!          do mj=1,m_def_fam
          do js=1,m_def_fam
      mj = def_family(js)
!            radnodj = nodes(mj)%horizon_size/horizon_factor
            radnodj = radnod(mj)
!>>02272009_YounDoh
            x1j = nodes(mj)%pos(1)
            x2j = nodes(mj)%pos(2)
            y1j = nodes(mj)%pos(1)+nodes(mj)%disp(1)
            y2j = nodes(mj)%pos(2)+nodes(mj)%disp(2)
            v1j = nodes(mj)%veloc(1)
            v2j = nodes(mj)%veloc(2)
            dvolj = nodes(mj)%volume
            nodtj = node_type(mj)
            xmassj = denst*dvolj

! >>031223

            ! Set the interaction distance for this pair of nodes (mi and mj).
            ! This is the threshold distance beyond which there is no unstructured interaction.

            ! Dist_init is the distance between the 2 nodes in the reference configuration.
            ! Dist_fac_init is the coefficient that multiplies the above to get the shortrange interaction
            ! distance for nodes initially far apart. This is usually > 1.

            ! Dist_nom is the nominal contact distance based on the max of the node radii.
            ! (The max is used so small nodes cannot sneak between large nodes.)
            ! Dist_fac_nom is the coefficient that multiplies the above. This is usually < 1.
            ! This determines the shortrange interaction distance for nodes initially close together.

            dist_init = dsqrt( (x1j-x1i)**2+(x2j-x2i)**2 )
!            dist_fac_init = max( shortrange_dist_fac_init(nodti), shortrange_dist_fac_init(nodtj) )
            dist_fac_init = 0.9d0
            dist_nom = 2.d0*dmax1(radnodi, radnodj)
!            dist_fac_nom = max( shortrange_dist_fac_nom(nodti), shortrange_dist_fac_nom(nodtj) )
            dist_fac_nom = 1.35d0
!            dist_fac_nom = 1.75d0

            deul = dmin1(dist_fac_nom*dist_nom, dist_fac_init*dist_init)

            !deul = 0.9*2.0*max(radnodi,radnodj)
            !deul = 1.5*deul
            !deul = min(deul, 0.9*dist_ref)

            ! Get the shortrange force factor (multiplies the peridynamic spring constant).
!            amp = 0.5*( shortrange_force_fac(nodti) + shortrange_force_fac(nodtj) )
!            amp = 15.0
            amp = 30.d0
! <<


            ! Find the current distance between this pair of nodes (mi and mj).
            distsq = (y1j-y1i)**2+(y2j-y2i)**2
            if(distsq<=deul**2.and.distsq>1.0d-6*deul**2) then

              ! Find the effective spring constant.
              ! Find the force on node mi due to unstructured interaction with node mj.
              fcof = One
              xmassij = 0.5d0*(xmassi+xmassj)
              call shortrange_force(y1i,y2i,y1j,y2j,deul,fcof,amp,v1i,v2i,v1j,v2j,xmassij,dvoli,dvolj,df1i,df2i,str,spcoef,dw)
              wt(mi) = wt(mi)+0.5d0*dw
              spsum(mi) = spsum(mi)+spcoef/denst
              nodes(mi)%new_veloc(1) = nodes(mi)%new_veloc(1) + df1i*dt/xmassi
              nodes(mi)%new_veloc(2) = nodes(mi)%new_veloc(2) + df2i*dt/xmassi


            end if
          end do
      end do

      return

end subroutine shortrange_motion
 
 
!-- shortrange_force
subroutine shortrange_force(y1a,y2a,y1b,y2b,deul,fcof,amp,v1a,v2a,v1b,v2b,xmass,dvola,dvolb,df1a,df2a,str,spcoef,dw)
 
      ! returns force on point a due to structureless interaction
      ! with point b.
 
      ! p is the scalar distance from a to b.
      p1 = y1b-y1a
      p2 = y2b-y2a
      p = dsqrt(p1**2+p2**2)
 
      ! find dp/dt.
      if(p>zero) then
        pdot = ((v1b-v1a)*(y1b-y1a)+(v2b-v2a)*(y2b-y2a))/p
      else
        pdot = zero
      end if
 
      ! amp is the amplification factor relative to the
      ! structured interaction forces
      !amp = 15.

      ! find the critical damping factor.
      ccrit = 2.d0*dsqrt(amp*fcof*dvola*dvolb*xmass)
      cdamp = 0.05d0*ccrit
 
      ! tension
      if(p>deul) then
        df1a = zero
        df2a = zero
        str = zero
        spcoef = zero
        dw = zero
 
      ! compression
      else if(p>=1.0d-10) then
        felast = amp*fcof*(p-deul)*dvola*dvolb
        fdamp = cdamp*pdot
        df1a = (felast+fdamp)*p1/p
        df2a = (felast+fdamp)*p2/p
        str = p-deul
        spcoef = amp*fcof*dvola
        dw = 0.5d0*amp*fcof*(p-deul)**2*dvola*dvolb
      else
        df1a = zero
        df2a = zero
        str = zero
        spcoef = zero
        dw = zero
      end if
 
      return
 
end subroutine shortrange_force

!****************************************************************************
end module peridynamic_computations
!****************************************************************************
