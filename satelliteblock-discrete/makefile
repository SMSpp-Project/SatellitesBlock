##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SatelliteBlock                                               #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by SatelliteBlock, i.e., core SMS++. These are          #
#   *not* copied into $(SATINC): adding those -I directives to the compile   #
#   commands will have to done by whatever "main" makefile is using this.    #
#   Analogously, any external library and the corresponding -L< libdirs >    #
#   will have to be added to the final linking command by  whatever "main"   #
#   makefile is using this.                                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by SatelliteBlock.                                                  #
#                                                                            #
#   The makefile defines internally (cf. MCFClssSlvr below) which            #
#   MCFSolver< :MCFClass > will be available.                                #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(SATSDR)    = the directory where the source is                 #
#                                                                            #
#   Output: $(SATOBJ)    = the final object(s) / library                   	 #
#           $(SATH)      = the .h files to include                      	 #
#           $(SATINC)    = the -I$( source directory )                	     #
#                                                                            #
#                              Luca Mencarelli	                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

SATOBJ = $(SATSDR)/obj/SatelliteBlock.o \
	$(SATSDR)/obj/ConstellationBlock.o \
	$(SATSDR)/obj/SatelliteSolver.o

SATINC = -I$(SATSDR)/include

SATH   = $(SATSDR)/include/SatelliteBlock.h \
	$(SATSDR)/include/ConstellationBlock.h \
	$(SATSDR)/include/SatelliteSolver.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SATOBJ) $(SATSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

all: $(SATOBJ)
#all : $(SATSDR)/obj/SatelliteBlock.o

$(SATSDR)/obj/SatelliteBlock.o: $(SATSDR)/src/SatelliteBlock.cpp \
	$(SATSDR)/include/SatelliteBlock.h $(SMS++OBJ)
	$(CC) -c $(SATSDR)/src/SatelliteBlock.cpp -o $@ \
	$(SATINC) $(SMS++INC) $(SW)

$(SATSDR)/obj/ConstellationBlock.o: \
	$(SATSDR)/src/ConstellationBlock.cpp $(SATH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SATSDR)/src/ConstellationBlock.cpp -o $@ \
	$(SATINC) $(SMS++INC) -I$(MILPINC) $(SW)

$(SATSDR)/obj/SatelliteSolver.o: \
	$(SATSDR)/src/SatelliteSolver.cpp $(SATH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SATSDR)/src/SatelliteSolver.cpp -o $@ \
	$(SATINC) $(SMS++INC) -I$(MILPINC) $(SW)

########################## End of makefile ###################################
