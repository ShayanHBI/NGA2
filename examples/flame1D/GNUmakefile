# NGA location if not yet defined
NGA_HOME ?= ../..

# Compilation parameters
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_HYPRE = TRUE
USE_LAPACK= TRUE
USE_FFTW  = TRUE
PROFILE   = FALSE
DEBUG     = FALSE
COMP      = gnu
EXEBASE   = nga

# Directories that contain user-defined code
Udirs   := src

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# External libraries are defined in .profile/.bashrc/.zshrc, but could be defined here as well

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
# Bdirs   := core grid config libraries data transform solver variable_density dvode finite_chem
Bdirs   := core finite_chem dvode variable_density data transform solver config grid libraries
Bpack   += $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
include $(Bpack)

# Inform user of Make.packages used
ifdef Ulocs
   $(info Taking user code from: $(Ulocs))
endif
$(info Taking base code from: $(Bdirs))

# Ignore MPI send/recv argument type mismatch error
F90FLAGS += -fallow-argument-mismatch

# Target definition
all: $(executable)
	@echo COMPILATION SUCCESSFUL

# NGA compilation rules
include $(NGA_HOME)/tools/GNUMake/Make.rules
