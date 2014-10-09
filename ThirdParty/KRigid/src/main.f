        PROGRAM main

!***********************************************************************
!        Written by Kumar Mithraratne
!        © Auckland Bioengineering Institute
!        V1.1 July 2013
!***********************************************************************
        
        IMPLICIT NONE
        
        INCLUDE 'rigidbody1.cmn'
        
        REAL*8 density(nbm),rog(nbm,3),volume(nbm)
        REAL*8 posni(nbm,3),lveli(nbm,3),lacci(nbm,3)
        REAL*8 ornti(nbm,3,3),aveli(nbm,3),aacci(nbm,3)
        REAL*8 eforce(nbm,3),etorque(nbm,3),bforce(nbm,3)
        REAL*8 time_initial,time_final,dt,output_interval
        REAL*8 translation(nbm,0:time_steps,9)
        REAL*8 rotation(nbm,0:time_steps,18)                
        !INTEGER 
        CHARACTER*100 rbname(nbm)
        
        print *, 'Starting'         
        CALL inertial_properties(rbname,volume,density,rog)

        print *, 'Init cond'         
        CALL initial_conditions(posni,lveli,lacci,ornti,aveli,aacci)
        
        print *, 'Extern force'
        CALL external_forces(eforce,etorque)
        
        print *, 'Body force'
        CALL body_forces(bforce)        
        
        print *, 'Time'
        CALL time_parameters(time_initial,time_final,dt,
     '       output_interval)
        
        print *, 'translate'
        CALL translational_motion(volume,density,
     '       posni,lveli,lacci,eforce,bforce,
     '       time_initial,time_final,dt,output_interval,translation)
        
        print *, 'rotate'
        CALL rotational_motion(volume,density,rog,
     '       ornti,aveli,aacci,etorque,
     '       time_initial,time_final,dt,output_interval,rotation)
                
        print *, 'output'
        CALL output(rbname,translation,rotation,output_interval)        
        
                
        STOP
        END
