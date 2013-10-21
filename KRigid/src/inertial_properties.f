        SUBROUTINE inertial_properties(rbname,volume,density,rog)        
        
        
!***********************************************************************
!        Written by Kumar Mithraratne
!        © Auckland Bioengineering Institute
!        V1.1 July 2013
!***********************************************************************
        
        IMPLICIT NONE
        
        INCLUDE 'rigidbody1.cmn'
        
        
        REAL*8 volume(nbm),density(nbm),rog(nbm,3)
        !INTEGER
        CHARACTER*100 rbname(nbm)
        
        !REAL*8 
        INTEGER nb,nj
        
                
        ! initialise
        DO nb=1,nbm,1
           !rbname(nb)=' '
           volume(nb)=0.0d0
           density(nb)=0.0d0
           DO nj=1,3,1
              rog(nb,nj)=0.0d0
           ENDDO
        ENDDO   
                
        !!!!!!!This data must be read in from an input file
        !shape - elliptical cylinder
        !length (L) - 150.0 mm
        !major axis (a) - 35.0 mm
        !minor axis (b) - 27.5 mm        
        !volume of elliptical cylinder (V) = pi*a*b*L - 453567.5 mm3
        !density of wood (rho)= 500 kg/m3 (average)
        !mass of elliptical cylinder (m) - 226.8 g
        !inertia (Ixx) = m*(0.25*b2+0.33*L2) - 1726869.3 gmm2 
        !nertia  (Iyy) = m*(0.25*a2+0.33*L2) = 1753447.5 gmm2
        !inertia (Izz) = 0.25*m*(a2+b2) = 112336.875 gmm2
        !RoGxx = sqrt((3*b2+4*L2)/12)
        !RoGyy = sqrt((3*a2+4*L2)/12)
        !RoGzz = sqrt((a2+b2)/4)
        
        nbt=2                !number of bodies (note: nbt in rigidbody1.cmn is set here)
        IF (nbt.GT.nbm) THEN
           WRITE(6,'(" Increase nbm in rigidbody1.cmn to",i6," !!!")')nbt
           STOP
        ENDIF
        
        DO nb=1,nbt,1
           rbname(nb)='projectile'
           volume(nb)=453567.5d0/cmm2cm   !volume m3
           density(nb)=500.0d0            !density of wood in kg/m3
           rog(nb,1)=87.687d0/smm2sm      !radius of gyration in m2
           rog(nb,2)=88.353d0/smm2sm 
           rog(nb,3)=22.256d0/smm2sm         
        ENDDO
        
             
        RETURN
        END
