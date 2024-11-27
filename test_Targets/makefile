##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for ConstellationBlock)                                #
#                                                                            #
#                              Luca Mencarelli                               #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = TARGETBlock_test

# basic directory
DIR = .

# debug switches
#SW = -g -glldb -fno-inline -std=c++17 -ferror-limit=1
# debug switches with address sanitizer and extra pedantic warning
#SW = -g3 -glldb -fno-inline -std=c++17 -ferror-limit=1 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter
# production switches
#SW = -O3 -std=c++17 -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION 
SW = -O3 -std=c++17 -DNDEBUG -DCLANG_1200_0_32_27_PATCH -ferror-limit=1

# compiler
CC = clang++ #-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(DIR)/$(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# SATBlock
SATSDR = ..
include $(SATSDR)/makefile-c

#LAGRANGIANDR = ../../lagrangiandualsolver
#BUNDLEDR = ../../bundlesolver

LgDSLVSDR = $(LAGRANGIANDR)
LgDSLVODR = $(LAGRANGIANDR)/obj/
include $(LAGRANGIANDR)/makefile-c

libNDOSDR = $(BUNDLEDR)/NdoFiOracle
include $(libNDOSDR)/lib/makefile

BNDSLVSDR = $(BUNDLEDR)
BNDSLVODR = $(BUNDLEDR)
include $(BUNDLEDR)/makefile

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(SATOBJ) 

# libraries
MLIB =  $(SATLIB) 

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o  
	$(CC) -o $(DIR)/$(NAME) $(SW) 

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(SATINC) 

# includes
MH =    $(SATH)

# compile command

$(DIR)/test.o: $(DIR)/test.cpp $(MH) 
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
