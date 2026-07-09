##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for SatellitesBlock)                                   #
#                                                                            #
#                                Donato Meoli                                #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# module name
NAME = SatellitesBlock_test

# basic directory
DIR = .

# compiler
CC = clang++

# common flags
COMMON_SW = -std=c++17 -ferror-limit=1 -Wno-enum-compare -DCLANG_1200_0_32_27_PATCH

# debug switches
SW_DEBUG = -g3 -glldb -fno-inline $(COMMON_SW)

# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# default target - build in release mode
default: SW = $(SW_RELEASE)
default: $(DIR)/$(NAME)

# debug target - build in debug mode
debug: SW = $(SW_DEBUG)
debug: $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# SMS++
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-inc

# SatellitesBlock
SBkSDR = ..
include $(SBkSDR)/makefile-s

LAGRANGIANDR = ../../LagrangianDualSolver
BUNDLEDR = ../../BundleSolver

libNDOSDR = $(BUNDLEDR)/NdoFiOracle
include $(BUNDLEDR)/NdoFiOracle/lib/makefile

BNDSLVSDR = $(BUNDLEDR)
BNDSLVODR = $(BUNDLEDR)
include $(BUNDLEDR)/makefile

libOSIINC = -I/usr/local/include/coin
libOSILIB = -rpath /usr/local/lib/ -lOsiClp -lClp -lOsi -lCoinUtils -lOsiCpx

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(SBkOBJ) $(SMS++OBJ)

# libraries
MLIB =  $(SBkLIB) $(SMS++LIB) -L"/Library/gurobi1301/macos_universal2/lib" -lgurobi130

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o
	$(CC) -o $(DIR)/$(NAME) $(DIR)/test.o ../../MILPSolver/obj/MILPSolver.o ../../MILPSolver/obj/GRBMILPSolver.o ../../LagrangianDualSolver/obj/*.o ../../BundleSolver/obj/*.o $(libNDOOBJ) $(MOBJ) $(libOSILIB) $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(SMS++INC) $(SBkINC)

# includes
MH =    $(SMS++H) $(SBkH)

# compile command
$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $(DIR)/test.cpp -o $@ -I../../MILPSolver/include/ -I../../LagrangianDualSolver/include/ -I../../BundleSolver/include/ -I$(libNDOSDR)/OPTUtils -I../../BundleSolver/NdoFiOracle/NDOSlver $(libOSIINC) $(MINC) $(SW)

# clean target
clean::
	rm -f $(DIR)/*.o $(DIR)/*~
	rm -f $(DIR)/$(NAME)

# phony targets
.PHONY: default debug clean

############################ End of makefile #################################
