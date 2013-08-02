# Variables to control Makefile operation

FORTRAN = f77
FLAG1 = -c

rigidbody: main.o \
	   inertial_properties.o \
	   initial_conditions.o \
	   external_forces.o \
	   body_forces.o \
	   time_parameters.o \
	   translational_motion.o \
	   rotational_motion.o \
	   atimesb.o \
	   output.o
	   $(FORTRAN) -o rigidbody *.o -L. -llapack -lblas
		
main.o: main.f
	$(FORTRAN) -c main.f
inertial_properties.o: inertial_properties.f	
	$(FORTRAN) -c inertial_properties.f
initial_conditions.o: initial_conditions.f		
	$(FORTRAN) -c initial_conditions.f
external_forces.o: external_forces.f		
	$(FORTRAN) -c external_forces.f	
body_forces.o: body_forces.f		
	$(FORTRAN) -c body_forces.f			
time_parameters.o: time_parameters.f		
	$(FORTRAN) -c time_parameters.f	
translational_motion.o: translational_motion.f		
	$(FORTRAN) -c translational_motion.f	
rotational_motion.o: rotational_motion.f		
	$(FORTRAN) -c rotational_motion.f	
atimesb.o: atimesb.f		
	$(FORTRAN) -c atimesb.f	
output.o: output.f		
	$(FORTRAN) -c output.f	

clean:
	@echo " Removing all .o files and executables..."
	rm -f *.o *.exdata *.txt *.com rigidbody
