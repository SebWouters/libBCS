###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Setting up the variables
# -----------------------------------------------------------------------------
LIBNAME = libBCS.so
CPPSRC	= source/BCS.cpp\
           source/BCSCoreSolver.cpp\
           source/BCS_cfunctions.cpp

OBJ	= $(CPPSRC:.cpp=.o)

LIBS	= -larpack -liomp5 -lpthread

CC	= icc
CXX	= icpc

CFLAGS	= -fPIC -openmp #fPIC also OK for icpc
LDFLAGS	= -shared


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(LIBNAME)...'
	@echo
	$(MAKE) $(LIBNAME)
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(LIBNAME) has been built successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(LIBNAME):	makefile $(OBJ)
	@echo; echo "Linker: creating $(LIBNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -Wl,-soname,$(LIBNAME) -o $(LIBNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo -n $(LIBNAME)
	@rm -f $(LIBNAME)
	@echo 'Done.'

#-----------------------------------------------------------------------------
# Make the documentation
#----------------------------------------------------------------------------
doc:
	@doxygen .doc-config

# ====================== End of file 'makefile.in' ========================== #
