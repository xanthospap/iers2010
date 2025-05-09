# Compiler and flags
FC = gfortran
FFLAGS = -O2

# Library name
LIBNAME = libtesti10for.a

# Source files for the library (adjust names here)
LIB_SRC = ORTHO_EOP.F CNMTX.F FUNDARG.F RG_ZONT2.F UTLIBR.F PMSDNUT2.F
LIB_OBJ = $(LIB_SRC:.F=.o)

# Main program source files (adjust names here)
BIN_SRC = TEST_ORTHO_EOP.F TEST_RG_ZONT2.F TEST_UTLIBR.F
BIN_OBJ = $(BIN_SRC:.F=.o)
BINARIES = $(BIN_SRC:.F=)

# Default target
all: $(LIBNAME) $(BINARIES)

# Build the library from object files
$(LIBNAME): $(LIB_OBJ)
	ar rcs $@ $^

# Rule to build the binary explicitly
$(BINARIES): %: %.F $(LIBNAME)
	$(FC) $(FFLAGS) -o $@ $< -L. -ltesti10for

# Compile .F to .o (Fortran 77 with preprocessing)
%.o: %.F
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f *.o *.mod $(LIBNAME) $(BINARIES)

