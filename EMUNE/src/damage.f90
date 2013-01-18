!****************************************************************************
!****************************************************************************
!*                                      *
!*  ModuleFile : damage_computations                    *
!*  Purpose    : This module includes subroutines to calculate the      *
!*        damage of initial bonds structure              *
!*                                      *
!*  PROGRAMMER : Dr. Youn Doh Ha                      *
!*                                      *
!****************************************************************************
module damage_computations

  use Global_variables
  use Objects
  use Input_subroutines
  use bond_family_computations
  use precision

  implicit none

contains

subroutine initial_broken_bond
 
  real(8):: x1,x2,y1,y2,y3, f1b, f2b
  integer(4) :: mi, mj, nc, js
  type(line)::line1

   
  ! Loop over all the nodes.
  do mi=1,nnodes

    nbroke = 0

    x1=nodes(mi)%pos(1)
    y1=nodes(mi)%pos(2)
    delta=nodes(mi)%horizon_size
    ! Get the family of node mi (all the nodes within its horizon, delta).
    call get_family(mi)

    ! Loop over nodes in the family of node mi.
    do js=1,mfam
      mj = family(js)

      x2=nodes(mj)%pos(1)
      y2=nodes(mj)%pos(2)

      line1%x1 = initial_crack_x0
      line1%y1 = initial_crack_y0
      line1%x2 = initial_crack_x1
      line1%y2 = initial_crack_y1
          
      if(LineIntersection(x1,y1,x2,y2,line1) .eqv. .TRUE.)then
!        if (x2.ne.x1) then
!          y3 = (y2-y1)/(x2-x1)*(line1%x1-x1)+y1
!          if (y3.le.line1%y2) then
!!          if(y1.le.line1%y2.or.y2.le.line1%y2) then
            broke(js,mi) = 1
            nbroke = nbroke + 1
!          endif
!        endif
      endif
    enddo

    f1b=dble(nbroke)
    f2b=dble(mfam)
    nodes(mi)%damage_index = f1b/f2b ! IF index is 0 then the node has no broken bond
    damage_index(mi)=f1b/f2b

  enddo

!  do mi=1,nnodes
!    nc = 0
!    do mj=1,nnodes
!      if (broke(mi,mj).eqv. .TRUE.1) nc = nc + 1
!    enddo
!    if (nc.gt.0) write(666,*) mi, nc
!  enddo

end subroutine initial_broken_bond
 
subroutine compute_critical_strain
 
  real(8):: h, syoung
  integer(4):: mi

  syoung = dsqrt(young)

  ! Loop over all the nodes.
  do mi=1,nnodes

    h = nodes(mi)%horizon_size
    if (sw_micro == 0) then
!    For constant micro-modulus
!      critical_strain(mi) = sqrt(8.d0*pi*fracture_energy/27.d0/young/nodes(mi)%horizon_size)
!      critical_strain(mi) = dsqrt(8.d0*pi*fracture_energy/27.d0/h)
      critical_strain(mi) = 2.d0/3.d0*dsqrt(2.d0*pi*fracture_energy/3.d0/h)
      critical_strain(mi) = critical_strain(mi)/syoung
    elseif (sw_micro == 1) then
!    For conical micro-modulus
!      critical_strain(mi) = dsqrt(10.d0*pi*fracture_energy/27.d0/h)
      critical_strain(mi) = dsqrt(10.d0*pi*fracture_energy/3.d0/h)/3.d0
      critical_strain(mi) = critical_strain(mi)/syoung
    endif

  enddo

end subroutine compute_critical_strain

real function CrossProduct(x1,y1,x2,y2)

  !Input: (x1,y1) and (x2,y2): (2D Vectors)
  real(8):: x1,y1,x2,y2
  real(8) :: value
  
  CrossProduct = (x1*y2 - x2*y1)

end function CrossProduct

!********************************************************************


logical function LineIntersection(Ax,Ay,Bx,By,Line1)

  !Input:
  !          C *
  !            |  
  !            |
  !        A *---|--------* B
  !            |  
  !            |
  !            * D      

  !Output: True or false

  real(8)::Ax,Ay,Bx,By,Cx,Cy,Dx,Dy
  type(line)::Line1
  real(8):: ABx,ABy,ACx,ACy,ADx,ADy,CDx,CDy,CBx,CBy,CAx,CAy
  real(8):: value1,value2,value3,value4

  Cx  = Line1%x1
  Cy  = Line1%y1
  Dx  = Line1%x2
  Dy  = Line1%y2

    ABx = Bx- Ax
  ABy = By- Ay

  ACx = Cx - Ax
  ACy = Cy - Ay

  ADx = Dx - Ax
  ADy = Dy - Ay

  CDx = Dx- Cx
  CDy = Dy- Cy

  CBx = Bx - Cx
  CBy = By - Cy

  CAx = Ax - Cx
  CAy = Ay - Cy

  !(ab x ac)*(ab x ad ) < 0  and (cd x ca)*(cd x cb) <0 
  value1 = CrossProduct(ACx,ACy,ABx,ABy)
  value2 = CrossProduct(ADx,ADy,ABx,ABy)
  value3 = CrossProduct(CAx,CAy,CDx,CDy)
  value4 = CrossProduct(CBx,CBy,CDx,CDy)

    LineIntersection = .false.
  if( (((value1) * (value2))<=zero) .and.  (((value3)*(value4))<=zero) )then
    LineIntersection = .true.
  endif
  




end function LineIntersection


!************************************************************************
!  If the body as an initial crack defined in line1
!  This subroutine will break all the bonds in the that cross this inital straight crack
!************************************************************************
subroutine broke_reg(mi, mj)
 
real(8):: x1,x2,y1,y2
integer(4) :: mi, mj
type(line)::line1

   
    x1=nodes(mi)%pos(1)
    y1=nodes(mi)%pos(2)
  x2=nodes(mj)%pos(1)
  y2=nodes(mj)%pos(2)
       
!    brokij=.FALSE. >>03062009_YounDoh

  line1%x1=initial_crack_x0
  line1%y1 = initial_crack_y0
  line1%x2= initial_crack_x1
  line1%y2 = initial_crack_y1
        
    if(LineIntersection(x1,y1,x2,y2,line1) .eqv. .TRUE.)then
    brokij=.TRUE.
  endif

end subroutine broke_reg

!****************************************************************************
end module damage_computations
!****************************************************************************
