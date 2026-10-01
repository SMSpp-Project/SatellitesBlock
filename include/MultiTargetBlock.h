/*--------------------------------------------------------------------------*/
/*------------------------ File MultiTargetBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiTargetBlock, which implements
 * the Block concept [see Block.h] for the Satellite Constellation Design
 * Problem decomposed by target, as a set of SingleTargetBlock linked by
 * consistency constraints.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MultiTargetBlock
 #define __MultiTargetBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "SingleTargetBlock.h"

#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*--------------------- MultiTargetBlock-RELATED TYPES ---------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlock_TYPES MultiTargetBlock-related types
 *
 * "Import" basic types from SingleTargetBlock.
 *  @{ */

using CNumber = SingleTargetBlock::CNumber;   ///< type of the costs
using c_RHSValue = RowConstraint::c_RHSValue; ///< type of the RHS values
using Vec_CNumber = SingleTargetBlock::Vec_CNumber; ///< vector of CNumber
using FNumber = SingleTargetBlock::FNumber;         ///< type of the values
using Vec_FNumber = SingleTargetBlock::Vec_FNumber; ///< vector of FNumber

using FMultiVector = std::vector< Vec_FNumber >;  ///< vector of Vec_FNumber
using CMultiVector = std::vector< Vec_CNumber >;  ///< vector of Vec_CNumber
using MultiSubset = std::vector< Block::Subset >; ///< vector of Subset

using Vec_Bool = std::vector< bool >; ///< vector of bool

/** @} end( group( MultiTargetBlock_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlock_CLASSES Classes in MultiTargetBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MultiTargetBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Satellite Constellation Design Problem, decomposed by target
/** MultiTargetBlock decomposes the constellation design problem *by
 * target* rather than by satellite (as ConstellationBlock does): it is
 * composed of one SingleTargetBlock [see SingleTargetBlock.h] per target,
 * each holding its own full copy of the satellite orbit-selection
 * Variable activation (there called pi[ j ][ c ], j indexing satellites
 * and c the candidate orbital configuration, common to all targets since
 * the candidate orbits [C] are computed once for the whole constellation,
 * see load()) and observability threshold theta[ j ]. Since physically
 * there is only one orbital configuration and one threshold per satellite,
 * shared by all targets, MultiTargetBlock links every pair of consecutive
 * SingleTargetBlock with "duplicate" consistency constraints
 *
 * \f[
 *  \pi_i[ j ][ c ] = \pi_{i+1}[ j ][ c ] , \quad \forall i \in
 *  \{ 1 , ... , m - 1 \} , \; j \in [s] , \; c \in [C]
 * \f]
 * \f[
 *  \theta_i[ j ] = \theta_{i+1}[ j ] , \quad \forall i \in
 *  \{ 1 , ... , m - 1 \} , \; j \in [s]
 * \f]
 *
 * where m is the number of targets and [s] = \{ 1 , 2 , ... , s \} is the
 * set of the satellites. Since
 * MultiTargetBlock defines no Objective of its own, and each
 * SingleTargetBlock's own objective is (a normalized fraction of) its
 * worst-case target revisit time, the overall problem minimizes the
 * average worst-case revisit time over all targets, subject to every
 * SingleTargetBlock agreeing on the same satellite configuration. */

class MultiTargetBlock : public Block {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of MultiTargetBlock
 /** Constructor of MultiTargetBlock. It accepts a pointer to the father
  * Block, which can be of any type, defaulting to nullptr so that this can
  * also be used as the void constructor. */

 explicit MultiTargetBlock( Block * father = nullptr )
  : Block( father ) , AR( false ) , satellites( 0 ) , targets( 0 ) ,
    time_step( 0 ) , horizon( 0 ) , indexOrbit( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: deletes the sub-Block and the abstract representation

 virtual ~MultiTargetBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance from the file with the given name
 /** Opens the file with the given name and loads the instance out of it
  * with load( std::istream & ); exception is thrown if the file cannot be
  * opened. */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// loads the instance out of an istream
 /** Loads the instance out of an istream, which has the format (comments
  * allowed, see eatcomments)
  *
  *     < time horizon, in hours >
  *     < length of the time step, in seconds >
  *     < number m of targets >
  *     m lines < latitude > < longitude >
  *     < number of satellites >
  *
  * with latitude and longitude of the targets in degrees. Then the
  * candidate orbits, the same for all the satellites, and the distances of
  * their ground tracks from the targets are computed, and one
  * SingleTargetBlock per target is created. The orbits are those of the
  * lowest altitude, in the range [ 400 , 1400 ] km, of the circular orbits
  * making an integer number of revolutions within the time horizon. Any
  * previous instance, including the sub-Block, is deleted. Exception is
  * thrown if the input is not correct or if there is no such altitude. If
  * there is any Solver attached to this MultiTargetBlock then a
  * NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the MultiTargetBlock
 /** The MultiTargetBlock has no Variable of its own: this just generates
  * those of all its SingleTargetBlock. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the MultiTargetBlock
 /** Generates the Variable and the Constraint of all the SingleTargetBlock
  * and then the "duplicate" constraints of the class comments, linking
  * them together. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR PRINTING & SAVING THE MultiTargetBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the MultiTargetBlock
 *  @{ */

 /// print the MultiTargetBlock on an ostream
 /** Prints the size of the instance (time horizon, time step, targets,
  * satellites, candidate orbits); vlvl is ignored. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/** @} ---------------------------------------------------------------------*/
/*---------- Methods for reading the data of the MultiTargetBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the MultiTargetBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns the time horizon [s]

 Index get_horizon( void ) const { return( horizon ); }

/*--------------------------------------------------------------------------*/
 /// returns the time step [s]

 Index get_timeStep( void ) const { return( time_step ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of satellites

 Index get_numSat( void ) const { return( satellites ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of targets

 Index get_numTarget( void ) const { return( targets ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of time stamps

 Index get_numTime( void ) const { return( horizon / time_step ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of candidate orbits

 Index get_numOrbits( void ) const { return( indexOrbit ); }

/*--------------------------------------------------------------------------*/
 /// returns the value of the Variable of orbit j of satellite i in the
 /// SingleTargetBlock of target k

 double get_activations( Index k , Index i , Index j ) const {
  return( static_cast< SingleTargetBlock * >( v_Block[ k ] )
            ->get_activation( i , j ) );
  }

/*--------------------------------------------------------------------------*/
 /// sets and fixes the Variable of orbit j of satellite i in the
 /// SingleTargetBlock of target k

 void set_activations( Index k , Index i , Index j , int value ) {
  static_cast< SingleTargetBlock * >( v_Block[ k ] )
   ->set_activation( i , j , value );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 bool AR; ///< true if the "duplicate" constraints have been generated

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/
/*----------------------------- PRIVATE FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h; // insert it in the Block factory

/*--------------------------------------------------------------------------*/

 Index satellites;  ///< the number of satellites
 Index targets;     ///< the number of targets
 FNumber time_step; ///< the time discretization step [s]
 FNumber horizon;   ///< the time horizon [s]
 Index indexOrbit;  ///< the number of candidate orbits

 boost::multi_array< FRowConstraint , 3 > duplicate_pi;
 ///< duplicate_pi[ i ][ j ][ k ], see the class comments
 boost::multi_array< FRowConstraint , 2 > duplicate_theta;
 ///< duplicate_theta[ i ][ j ], see the class comments

/*--------------------------------------------------------------------------*/

 }; // end( class( MultiTargetBlock ) )

/** @} end( group( MultiTargetBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* MultiTargetBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiTargetBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
