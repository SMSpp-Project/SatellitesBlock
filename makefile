##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SatellitesBlock                                              #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by            #
#   SatellitesBlock.                                                         #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(SMS++INC)  = the -I$( core SMS++ directory )                   #
#           $(SMS++OBJ)  = the libSMS++ library itself                       #
#               $(SBkSDR)  = the directory where the source is               #
#                                                                            #
#   Output: $(SBkOBJ)  = the final object(s) / library                       #
#                    $(SBkH)    = the .h files to include                    #
#                  $(SBkINC)  = the -I$( source directory )                  #
#                                                                            #
#                             Antonio Frangioni                              #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

SBkOBJ = $(SBkSDR)/obj/ConstellationBlock.o \
         $(SBkSDR)/obj/MultiTargetBlock.o \
         $(SBkSDR)/obj/MultiTargetBlockv2.o \
         $(SBkSDR)/obj/SatelliteBlock.o \
         $(SBkSDR)/obj/SatelliteSolver.o \
         $(SBkSDR)/obj/SingleTargetBlock.o \
		 $(SBkSDR)/obj/ConstellationBlock_discrete.o \
		 $(SBkSDR)/obj/SatelliteBlock_discrete.o \
         $(SBkSDR)/obj/SatelliteSolver_discrete.o 

SBkINC = -I$(SBkSDR)/include

SBkH   = $(SBkSDR)/include/ConstellationBlock.h \
         $(SBkSDR)/include/MultiTargetBlock.h \
         $(SBkSDR)/include/MultiTargetBlockv2.h \
         $(SBkSDR)/include/SatelliteBlock.h \
         $(SBkSDR)/include/SatelliteSolver.h \
         $(SBkSDR)/include/SingleTargetBlock.h \
		 $(SBkSDR)/include/ConstellationBlock_discrete.h \
		 $(SBkSDR)/include/SatelliteBlock_discrete.h \
         $(SBkSDR)/include/SatelliteSolver_discrete.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SBkOBJ) $(SBkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(SBkSDR)/obj/ConstellationBlock.o: $(SBkSDR)/src/ConstellationBlock.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/ConstellationBlock.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/MultiTargetBlock.o: $(SBkSDR)/src/MultiTargetBlock.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/MultiTargetBlock.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/MultiTargetBlockv2.o: $(SBkSDR)/src/MultiTargetBlockv2.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/MultiTargetBlockv2.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/SatelliteBlock.o: $(SBkSDR)/src/SatelliteBlock.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/SatelliteBlock.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/SatelliteSolver.o: $(SBkSDR)/src/SatelliteSolver.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/SatelliteSolver.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/SingleTargetBlock.o: $(SBkSDR)/src/SingleTargetBlock.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/SingleTargetBlock.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/ConstellationBlock_discrete.o: $(SBkSDR)/src/ConstellationBlock_discrete.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/ConstellationBlock_discrete.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/SatelliteBlock_discrete.o: $(SBkSDR)/src/SatelliteBlock_discrete.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/SatelliteBlock_discrete.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

$(SBkSDR)/obj/SatelliteSolver_discrete.o: $(SBkSDR)/src/SatelliteSolver_discrete.cpp $(SBkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SBkSDR)/src/SatelliteSolver_discrete.cpp -o $@ \
	$(SBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
