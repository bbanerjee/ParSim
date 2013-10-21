        SUBROUTINE time_parameters(time_initial,time_final,dt,
     '             output_interval)
        
        
!***********************************************************************
!        Written by Kumar Mithraratne
!        © Auckland Bioengineering Institute
!        V1.1 July 2013
!***********************************************************************
        
        IMPLICIT NONE
        
        INCLUDE 'rigidbody1.cmn'
        
        
        REAL*8 time_initial,time_final,dt,output_interval
        !INTEGER
        
        !REAL*8 
        !INTEGER 

         
        !!!!!!!This information must be read in from an input file
        
        time_initial=0.0d0;      !seconds
        time_final=1.0d0;        !seconds
        dt=0.0001d0;             !seconds
        output_interval=0.05d0   !seconds
        

        RETURN
        END


