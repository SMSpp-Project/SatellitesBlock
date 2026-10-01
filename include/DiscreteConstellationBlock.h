/*--------------------------------------------------------------------------*/
/*------------------- File DiscreteConstellationBlock.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteConstellationBlock, which
 * implements the Block concept [see Block.h] for the Satellite
 * Constellation Design Problem in the form where the observability
 * threshold of each satellite is discretized into a finite number of
 * levels.
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

#ifndef __DiscreteConstellationBlock
 #define __DiscreteConstellationBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "DiscreteSatelliteBlock.h"

#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup DiscreteConstellationBlock_CLASSES Discrete constellation
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS DiscreteConstellationBlock --------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Satellite Constellation Design Problem, discretized form
/** DiscreteConstellationBlock is the counterpart of ConstellationBlock [see
 * ConstellationBlock.h] in which the observability threshold
 * \f$ \theta^{\max} \f$ of each satellite, a continuous Variable in
 * SatelliteBlock, is discretized into ell = 3 levels, the same for all the
 * satellites: level l is \f$ \theta_l = ( ell - l ) / ell \cdot
 * \Theta_{\min} / 10 \f$, where \f$ \Theta_{\min} \f$ is the angular
 * resolution computed by load(), so that level 0 is the loosest. Its
 * sub-Block are one DiscreteSatelliteBlock [see DiscreteSatelliteBlock.h]
 * per satellite.
 *
 * load() reads the instance, builds the candidate orbits of each satellite
 * with the same ground-track propagation as ConstellationBlock::load(), and
 * computes once and for all the 0/1 coefficient
 *
 *     obs[ k ][ j ][ i ][ o ][ l ]
 *
 * which is 1 if satellite k, placed on its candidate orbit o with the
 * threshold level l, observes target i at time stamp j, i.e., if both the
 * latitude and the longitude distance between the ground track and the
 * target are within the threshold \f$ \theta_l \f$. The observability
 * conditions, that SatelliteBlock expresses by big-M constraints, are then
 * linear constraints in the Variable y[ o ][ l ] of the
 * DiscreteSatelliteBlock, with the obs as coefficients. With [s] the set of
 * the satellites and pp = T / ( dt periods[ i ] ) the number of time stamps
 * of a revisit window of target i, these constraints are
 *
 * \f[
 *  \sum_{ k \in [s] } \sum_{ o , l } \Bigl( \sum_{ j' = j\,pp }^{
 *  (j+1)pp - 1 } obs[ k ][ j' ][ i ][ o ][ l ] \Bigr) y_k[ o ][ l ] \geq 1
 *  \quad j = 0 , \ldots , periods[ i ] - 1                         \qquad (1)
 * \f]
 *
 * i.e., every target is observed at least once in each of its revisit
 * windows (the "observation" constraints),
 *
 * \f[
 *  \sum_{ k \in [s] } \sum_{ o , l } obs[ k ][ j ][ i ][ o ][ l ]
 *  y_k[ o ][ l ] \leq 1                                            \qquad (2)
 * \f]
 *
 * for every target i and time stamp j, i.e., no two satellites observe the
 * same target at the same time (the "observation1" constraints), and
 *
 * \f[
 *  \sum_{ k \in [s] } \sum_{ o , l } ( \theta_l - 0.9 \theta_0 )
 *  y_k[ o ][ l ] \leq 0                                            \qquad (3)
 * \f]
 *
 * i.e., the average threshold of the active satellites is at most 0.9
 * times the loosest level \f$ \theta_0 \f$ (the "thetaM" constraint). The
 * DiscreteConstellationBlock has no Objective of its own: the objective of
 * the problem is the sum of those of the DiscreteSatelliteBlock, i.e., the
 * number of active satellites. */

class DiscreteConstellationBlock : public Block {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * "Import" the basic types from DiscreteSatelliteBlock.
 *  @{ */

 using FNumber = DiscreteSatelliteBlock::FNumber; ///< type of the values
 using CNumber = DiscreteSatelliteBlock::CNumber; ///< type of the costs
 using Vec_FNumber = DiscreteSatelliteBlock::Vec_FNumber;
 ///< a vector of FNumber
 using Vec_CNumber = DiscreteSatelliteBlock::Vec_CNumber;
 ///< a vector of CNumber

/** @} ---------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of DiscreteConstellationBlock
 /** Constructor of DiscreteConstellationBlock. It accepts a pointer to the
  * father Block, which can be of any type, defaulting to nullptr so that
  * this can also be used as the void constructor. */

 explicit DiscreteConstellationBlock( Block * father = nullptr )
  : Block( father ) , satellites( 0 ) , targets( 0 ) , time_step( 0 ) ,
    horizon( 0 ) , AR( false ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: deletes the sub-Block and the abstract representation

 virtual ~DiscreteConstellationBlock() { guts_of_destructor(); }

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
 /** Loads the instance out of an istream, which has the format of the input
  * of ConstellationBlock::load() (comments allowed, see eatcomments):
  *
  *     < time horizon, in hours >
  *     < length of the time step, in seconds >
  *     < number m of targets >
  *     m lines < number of revisit periods > < latitude > < longitude >
  *     < number of satellites >
  *
  * with latitude and longitude of the targets in degrees. Then the
  * candidate orbits and the obs[][][][][] coefficients are computed (see
  * the class comments), and one DiscreteSatelliteBlock per satellite is
  * created; the satellites are split in three groups of (almost) the same
  * size, each one using a different altitude among those, in the range
  * [ 400 , 1400 ] km, of the circular orbits making an integer number of
  * revolutions within the time horizon. Any previous instance, including
  * the sub-Block, is deleted. Exception is thrown if the input is not
  * correct or if the time horizon allows less than three altitudes. If
  * there is any Solver attached to this DiscreteConstellationBlock then a
  * NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the DiscreteConstellationBlock
 /** The DiscreteConstellationBlock has no Variable of its own: this just
  * generates those of all its DiscreteSatelliteBlock. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the DiscreteConstellationBlock
 /** Generates the Variable and the Constraint of all the
  * DiscreteSatelliteBlock and then the constraints (1), (2) and (3) of the
  * class comments, linking them together. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*---------- Methods for reading the data of the DiscreteConstellationBlock */
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the DiscreteConstellationBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of satellites

 Index get_numSat( void ) const { return( satellites ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of candidate orbits of satellite k

 Index get_orbits( Index k ) const {
  return(
   static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->get_orbits() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of threshold levels of satellite k

 Index get_ell( Index k ) const {
  return(
   static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->get_ell() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable y[ o ][ l ] of satellite k

 FNumber get_solution( Index k , Index o , Index l ) const {
  return(
   static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->get_y( o , l ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution satisfies the linking constraints
 /** Returns true if the solution encoded in the current value of the
  * Variable of the DiscreteSatelliteBlock satisfies the constraints (1),
  * (2) and (3) of the class comments, which requires that they have been
  * generated (false is returned otherwise); the constraints of each
  * DiscreteSatelliteBlock are not checked, this being done by their own
  * is_feasible(). The tolerance, by default 1e-1 and relative, is taken
  * out of fsbc, or, if this is not a valid Configuration, out of
  * f_BlockConfig->f_is_feasible_Configuration; a valid Configuration is
  * either a SimpleConfiguration< double > containing the tolerance, or a
  * SimpleConfiguration< std::pair< double , int > > containing the
  * tolerance and whether it is relative (nonzero) or absolute (zero). */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------ METHODS FOR PRINTING & SAVING THE DiscreteConstellationBlock ------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the DiscreteConstellationBlock
 *  @{ */

 /// print the DiscreteConstellationBlock on an ostream
 /** Prints the size of the instance (time horizon, time step, targets,
  * satellites) and the number of candidate orbits and threshold levels of
  * each satellite; vlvl is ignored. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/** @} ---------------------------------------------------------------------*/
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

 std::vector< Index > periods; ///< the number of revisit periods per target
 std::vector< Index > indexOrbitSat; ///< the number of orbits per satellite
 std::vector< Index > ellSat;        ///< the number of levels per satellite

 boost::multi_array< double , 5 > obs;
 ///< obs[ k ][ j ][ i ][ o ][ l ], see the class comments
 boost::multi_array< double , 2 > thetaValSat;
 ///< thetaValSat[ k ][ l ], threshold l of satellite k

 FRowConstraint thetaM; ///< the constraint (3)

 boost::multi_array< FRowConstraint , 2 > observation;
 ///< the constraints (1)
 boost::multi_array< FRowConstraint , 2 > observation1;
 ///< the constraints (2)

 bool AR; ///< true if the constraints (1), (2) and (3) have been generated

/*--------------------------------------------------------------------------*/

 }; // end( class( DiscreteConstellationBlock ) )

/** @} end( group( DiscreteConstellationBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* DiscreteConstellationBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------- End File DiscreteConstellationBlock.h ------------------*/
/*--------------------------------------------------------------------------*/
