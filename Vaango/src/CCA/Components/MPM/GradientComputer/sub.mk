# Makefile fragment for this subdirectory 


SRCDIR   := CCA/Components/MPM/GradientComputer

SRCS     += \
	$(SRCDIR)/DeformationGradientComputer.cc \
	$(SRCDIR)/DisplacementGradientComputer.cc \
	$(SRCDIR)/GradientComputer.cc \
	$(SRCDIR)/VelocityGradientComputer.cc
