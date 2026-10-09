/*--------------------------------------------------------------------------*/
/*----------------------- File ConstellationBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class ConstellationBlock, which implements
 * the Block concept [see Block.h] for the Satellite Constellation Design
 * Problem, as a set of SatelliteBlock linked by the observability
 * constraints of the targets.
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

#ifndef __ConstellationBlock
 #define __ConstellationBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "SatelliteBlock.h"

#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------- ConstellationBlock-RELATED TYPES --------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ConstellationBlock_TYPES ConstellationBlock-related types
 *
 * "Import" basic types from SatelliteBlock.
 *  @{ */

using CNumber = SatelliteBlock::CNumber;         ///< type of the costs
using c_RHSValue = RowConstraint::c_RHSValue;    ///< type of the RHS values
using Vec_CNumber = SatelliteBlock::Vec_CNumber; ///< vector of CNumber
using FNumber = SatelliteBlock::FNumber;         ///< type of the values
using Vec_FNumber = SatelliteBlock::Vec_FNumber; ///< vector of FNumber

using FMultiVector = std::vector< Vec_FNumber >;  ///< vector of Vec_FNumber
using CMultiVector = std::vector< Vec_CNumber >;  ///< vector of Vec_CNumber
using MultiSubset = std::vector< Block::Subset >; ///< vector of Subset

using Vec_Bool = std::vector< bool >; ///< vector of bool

/** @} end( group( ConstellationBlock_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ConstellationBlock_CLASSES Classes in ConstellationBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ConstellationBlock ------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Satellite Constellation Design Problem
/** ConstellationBlock represents the Satellite Constellation Design
 * Problem: choosing the orbits of a constellation of at most s satellites,
 * among a finite set of candidate ones, so that every target m is observed
 * by some satellite at least once in each of its revisit windows, with the
 * fewest active satellites. Its sub-Block are one SatelliteBlock [see
 * SatelliteBlock.h] per satellite, linked by the observability constraints
 *
 * \f[
 *  \sum_{ i \in [s] } \sum_{ t \in T( k , \Delta t_m , dt ) }
 *  \xi_i[ m ][ t ] \geq 1 \quad k = 1 , \ldots ,
 *  \lfloor T / \Delta t_m \rfloor \; , \; m \in \mathcal{X}      \qquad (1)
 * \f]
 *
 * where \f$ \mathcal{X} \f$ is the set of the targets, T the time horizon,
 * dt the time step, \f$ \Delta t_m \f$ the revisit time of target m,
 * \f$ T( k , \Delta t_m , dt ) \f$ the set of time stamps of its k-th
 * revisit window, [s] = { 1 , ... , s } the set of the satellites, and
 * \f$ \xi_i[ m ][ t ] \f$ the Variable of SatelliteBlock i that is 1 if the
 * satellite observes target m at time stamp t. Besides (1), the
 * "observation" constraints, ConstellationBlock has the "observation1"
 * constraints
 *
 * \f[
 *  \sum_{ i \in [s] } \xi_i[ m ][ t ] \leq 1                     \qquad (2)
 * \f]
 *
 * for every target m and time stamp t (no two satellites observe the same
 * target at the same time), and the "thetaM" constraint
 *
 * \f[
 *  \sum_{ i \in [s] } \theta_i \leq 0.9 \, \bar\theta
 *  \sum_{ i \in [s] } \zeta_i                                    \qquad (3)
 * \f]
 *
 * bounding the average observability threshold \f$ \theta_i \f$ of the
 * active satellites (those with \f$ \zeta_i = 1 \f$) with the reference
 * threshold \f$ \bar\theta \f$ computed by load(). ConstellationBlock has
 * no Objective of its own: the objective of the problem is the sum of those
 * of the SatelliteBlock, i.e., the number of active satellites. */

class ConstellationBlock : public Block {

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

 /// constructor of ConstellationBlock
 /** Constructor of ConstellationBlock. It accepts a pointer to the father
  * Block, which can be of any type, defaulting to nullptr so that this can
  * also be used as the void constructor. */

 explicit ConstellationBlock( Block * father = nullptr )
  : Block( father ) , AR( false ) , satellites( 0 ) , targets( 0 ) ,
    time_step( 0 ) , horizon( 0 ) , thetaValF( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: deletes the sub-Block and the abstract representation

 virtual ~ConstellationBlock() { guts_of_destructor(); }

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
  *     m lines < number of revisit periods > < latitude > < longitude >
  *     < number of satellites >
  *
  * with latitude and longitude of the targets in degrees. Then the
  * candidate orbits of each satellite and the distances of their ground
  * tracks from the targets are computed, and one SatelliteBlock per
  * satellite is created; the satellites are split in three groups of
  * (almost) the same size, each one using a different altitude among
  * those, in the range [ 400 , 1400 ] km, of the circular orbits making an
  * integer number of revolutions within the time horizon. Any previous
  * instance, including the sub-Block, is deleted. Exception is thrown if
  * the input is not correct or if the time horizon allows less than three
  * altitudes. If there is any Solver attached to this ConstellationBlock
  * then a NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the ConstellationBlock
 /** The ConstellationBlock has no Variable of its own: this just generates
  * those of all its SatelliteBlock. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the ConstellationBlock
 /** Generates the Variable and the Constraint of all the SatelliteBlock and
  * then the constraints (1), (2) and (3) of the class comments, linking
  * them together. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution satisfies the linking constraints
 /** Returns true if the solution encoded in the current value of the
  * Variable of the SatelliteBlock satisfies the constraints (1), (2) and
  * (3) of the class comments, which requires that they have been generated
  * (false is returned otherwise); the constraints of each SatelliteBlock
  * are not checked, this being done by their own is_feasible(). The
  * tolerance, by default 1e-1 and relative, is taken out of fsbc, or, if
  * this is not a valid Configuration, out of
  * f_BlockConfig->f_is_feasible_Configuration; a valid Configuration is
  * either a SimpleConfiguration< double > containing the tolerance, or a
  * SimpleConfiguration< std::pair< double , int > > containing the
  * tolerance and whether it is relative (nonzero) or absolute (zero). */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE ConstellationBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the ConstellationBlock
 *  @{ */

 /// print the ConstellationBlock on an ostream
 /** Prints the size of the instance (time horizon, time step, targets,
  * satellites) and the number of candidate orbits of each satellite; vlvl
  * is ignored. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------- Methods for reading the data of the ConstellationBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the ConstellationBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the threshold Variable of satellite k

 double get_thetaVar( Index k ) const { return( SB( k )->get_thetaVar() ); }

/*--------------------------------------------------------------------------*/
 /// returns the altitude of the orbits of satellite k

 double get_alt( Index k ) const { return( SB( k )->get_alt() ); }

/*--------------------------------------------------------------------------*/
 /// returns the half-aperture of the sensor of satellite k

 double get_alpha( Index k ) const { return( SB( k )->get_alpha() ); }

/*--------------------------------------------------------------------------*/
 /// returns the maximum threshold of satellite k

 double get_theta( Index k ) const { return( SB( k )->get_theta() ); }

/*--------------------------------------------------------------------------*/
 /// returns the latitude distance of orbit i of satellite k from target n
 /// at time stamp t

 double get_Lat( Index k , Index i , Index n , Index t ) const {
  return( SB( k )->get_Delta_lat( i , n , t ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the longitude distance of orbit i of satellite k from target n
 /// at time stamp t

 double get_Long( Index k , Index i , Index n , Index t ) const {
  return( SB( k )->get_Delta_long( i , n , t ) );
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
 /// returns the number of revisit periods of target k

 Index get_period( Index k ) const { return( periods[ k ] ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of candidate orbits of satellite k

 Index get_numOrbits( Index k ) const { return( SB( k )->get_numOrbit() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the activation Variable of satellite k

 double get_zs( Index k ) const { return( SB( k )->get_zeta() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the observation Variable of satellite k
 /// for target n at time stamp t

 double get_xis( Index k , Index n , Index t ) const {
  return( SB( k )->get_xi( n , t ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable of orbit j of satellite k

 double get_activations( Index k , Index j ) const {
  return( SB( k )->get_activation( j ) );
  }

/*--------------------------------------------------------------------------*/
 /// sets and fixes the Variable of orbit j of satellite k to value

 void set_activations( Index k , Index j , int value ) {
  SB( k )->set_activation( j , value );
  }

/*--------------------------------------------------------------------------*/
 /// sets and fixes the activation Variable of satellite k to value

 void set_zs( Index k , int value ) { SB( k )->set_zeta( value ); }

/** @} ---------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// returns the k-th sub-Block, as a SatelliteBlock

 SatelliteBlock * SB( Index k ) const {
  return( static_cast< SatelliteBlock * >( v_Block[ k ] ) );
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 bool AR; ///< true if the constraints (1), (2) and (3) have been generated

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
 double thetaValF;  ///< the reference threshold of the satellites

 Vec_CNumber periods; ///< the number of revisit periods per target

 boost::multi_array< FRowConstraint , 2 > observation;
 ///< the constraints (1)
 boost::multi_array< FRowConstraint , 2 > observation1;
 ///< the constraints (2)
 std::vector< FRowConstraint > thetaM; ///< the constraint (3)

/*--------------------------------------------------------------------------*/

 }; // end( class( ConstellationBlock ) )

/** @} end( group( ConstellationBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ConstellationBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File ConstellationBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
