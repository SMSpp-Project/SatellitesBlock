/*--------------------------------------------------------------------------*/
/*--------------------- File DiscreteSatelliteBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteSatelliteBlock, which
 * implements the Block concept [see Block.h] for a single satellite of the
 * Satellite Constellation Design Problem in the discretized form used by
 * DiscreteConstellationBlock, and for the class DiscreteSatelliteSolution
 * holding a solution of it.
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

#ifndef __DiscreteSatelliteBlock
 #define __DiscreteSatelliteBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
class DiscreteSatelliteBlock; // forward declaration of DiscreteSatelliteBlock

class DiscreteSatelliteSolution; // forward declaration of the Solution

/*--------------------------------------------------------------------------*/
/*------------------ DiscreteSatelliteBlock-RELATED TYPES ------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup DiscreteSatelliteBlock_TYPES DiscreteSatelliteBlock types
 *  @{ */

using p_DiscreteSatelliteBlock = DiscreteSatelliteBlock *;
///< a pointer to DiscreteSatelliteBlock

using Vec_DiscreteSatelliteBlock = std::vector< p_DiscreteSatelliteBlock >;
///< a vector of pointers to DiscreteSatelliteBlock

using Vec_DiscreteSatelliteBlock_it = Vec_DiscreteSatelliteBlock::iterator;
///< iterator for a Vec_DiscreteSatelliteBlock

using c_Vec_DiscreteSatelliteBlock = const Vec_DiscreteSatelliteBlock;
///< a const vector of pointers to DiscreteSatelliteBlock

using c_Vec_DiscreteSatelliteBlock_it =
 c_Vec_DiscreteSatelliteBlock::iterator;
///< iterator for a c_Vec_DiscreteSatelliteBlock

/** @} end( group( DiscreteSatelliteBlock_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup DiscreteSatelliteBlock_CLASSES Discrete satellite
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS DiscreteSatelliteBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a single satellite of a DiscreteConstellationBlock
/** DiscreteSatelliteBlock is the counterpart of SatelliteBlock [see
 * SatelliteBlock.h] in the discretized form of the Satellite Constellation
 * Design Problem handled by DiscreteConstellationBlock [see
 * DiscreteConstellationBlock.h]. In SatelliteBlock the observability
 * threshold \f$ \theta^{\max} \f$ of the satellite is a continuous Variable,
 * linked to the observation Variable by big-M constraints; here, instead,
 * \f$ \theta^{\max} \f$ is discretized into ell levels, so that the choice of
 * the satellite is that of a pair (orbit , threshold level) among the
 * OrbitSet x ell possible ones, represented by the binary Variable
 *
 * \f[
 *  y[ o ][ l ] \in \{ 0 , 1 \} \quad o = 0 , \ldots , OrbitSet - 1 \; , \;
 *  l = 0 , \ldots , ell - 1
 * \f]
 *
 * with y[ o ][ l ] = 1 if the satellite is placed on the candidate orbit o
 * with the threshold level l. The only constraint is that at most one pair
 * is chosen
 *
 * \f[
 *  \sum_{ o = 0 }^{ OrbitSet - 1 } \sum_{ l = 0 }^{ ell - 1 } y[ o ][ l ]
 *  \leq 1 \; ,
 * \f]
 *
 * the satellite being unused (inactive) if all the y are zero. The
 * objective is the LinearFunction with unit coefficients over all the y,
 * i.e., 1 if the satellite is active and 0 otherwise, so that the sum of
 * the objectives of the DiscreteSatelliteBlock of a
 * DiscreteConstellationBlock is the number of active satellites.
 *
 * Which targets each pair (o , l) observes, and at which time stamps, is not
 * known to the DiscreteSatelliteBlock: it is precomputed by
 * DiscreteConstellationBlock, which uses it as the coefficients of the y in
 * its own observability constraints. As a consequence, once these are
 * dualized the DiscreteSatelliteBlock is a trivial problem, solved "by
 * inspection" by DiscreteSatelliteSolver [see DiscreteSatelliteSolver.h]. */

class DiscreteSatelliteBlock : public Block {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * DiscreteSatelliteBlock defines three main public types:
 *
 * - FNumber, the type of the values of the Variable and of the data;
 *
 * - CNumber, the type of the objective coefficients;
 *
 * - FONumber, the type of the objective function value.
 *
 * They are all double; they are kept as separate names for uniformity with
 * SatelliteBlock and with the other Block of the module.
 *  @{ */

 using FNumber = double;          ///< type of the values
 using c_FNumber = const FNumber; ///< a read-only FNumber

 using Vec_FNumber = std::vector< FNumber >; ///< a vector of FNumber
 using c_Vec_FNumber = const Vec_FNumber;    ///< a const vector of FNumber

 using Vec_FNumber_it = Vec_FNumber::iterator; ///< iterator in Vec_FNumber
 using c_Vec_FNumber_it = Vec_FNumber::const_iterator;
 ///< const iterator in Vec_FNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 using CNumber = double;          ///< type of the costs
 using c_CNumber = const CNumber; ///< a read-only CNumber

 using Vec_CNumber = std::vector< CNumber >; ///< a vector of CNumber
 using c_Vec_CNumber = const Vec_CNumber;    ///< a const vector of CNumber

 using Vec_CNumber_it = Vec_CNumber::iterator; ///< iterator in Vec_CNumber
 using c_Vec_CNumber_it = Vec_CNumber::const_iterator;
 ///< const iterator in Vec_CNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 using FONumber = double;           ///< type of the objective function value
 using c_FONumber = const FONumber; ///< a read-only FONumber

 using Vec_FONumber = std::vector< FONumber >; ///< a vector of FONumber
 using c_Vec_FONumber = const Vec_FONumber;    ///< a const vector of FONumber

/** @} ---------------------------------------------------------------------*/
/*-------------------------------- FRIENDS ---------------------------------*/
/*--------------------------------------------------------------------------*/

 friend DiscreteSatelliteSolution; ///< make DiscreteSatelliteSolution friend

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of DiscreteSatelliteBlock, taking a pointer to the father
 /** Constructor of DiscreteSatelliteBlock. It accepts a pointer to the
  * father Block, which can be of any type, defaulting to nullptr so that
  * this can also be used as the void constructor. */

 explicit DiscreteSatelliteBlock( Block * father = nullptr )
  : Block( father ) , AR1( 0 ) , AR2( 0 ) , AR3( 0 ) , OrbitSet( 0 ) ,
    ell( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: deletes the abstract representation, if any

 virtual ~DiscreteSatelliteBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance out of the size of the (orbit , level) grid
 /** Loads the instance, which is entirely described by the number
  * n_orbits of candidate orbits and the number n_levels of levels of the
  * observability threshold (both are computed by
  * DiscreteConstellationBlock::load()). Any previous abstract
  * representation is deleted. If there is any Solver attached to this
  * DiscreteSatelliteBlock then a NBModification (the "nuclear option") is
  * issued. */

 void load( Index n_orbits , Index n_levels );

/*--------------------------------------------------------------------------*/
 /// loads the instance out of an istream
 /** Loads the instance out of an istream, which must contain (possibly
  * preceded by comments, see eatcomments) the two numbers
  *
  *     < number of candidate orbits > < number of threshold levels >
  *
  * that are then passed to load( Index , Index ). The format parameter
  * frmt is ignored. Exception is thrown if the two numbers cannot be read. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// loads the instance from the file with the given name
 /** The version of Block, which opens the file and calls
  * load( std::istream & ). */

 using Block::load;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the DiscreteSatelliteBlock
 /** Generates the OrbitSet x ell binary Variable y[][]. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the DiscreteSatelliteBlock
 /** Generates the single constraint sum_{ o , l } y[ o ][ l ] <= 1. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the Objective of the DiscreteSatelliteBlock
 /** Generates the FRealObjective whose LinearFunction has unit coefficient
  * on every y[][], in "row-major" order, i.e., that of y[ o ][ l ] is the
  * coefficient with index o * ell + l. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------- Methods for reading the data of the DiscreteSatelliteBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the DiscreteSatelliteBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of candidate orbits

 [[nodiscard]] Index get_orbits( void ) const { return( OrbitSet ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of levels of the observability threshold

 [[nodiscard]] Index get_ell( void ) const { return( ell ); }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the Variable y[ o ][ l ]
 /** Returns a pointer to the Variable y[ o ][ l ]; this clearly requires
  * that generate_abstract_variables() has been called. */

 [[nodiscard]] ColVariable * i2p_y( Index o , Index l ) const {
  return( const_cast< ColVariable * >( &y[ o ][ l ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable y[ o ][ l ]

 [[nodiscard]] FNumber get_y( Index o , Index l ) const {
  return( y[ o ][ l ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// sets the current value of the Variable y[ o ][ l ]

 void set_y( Index o , Index l , FNumber value ) {
  y[ o ][ l ].set_value( value );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the y[][]
  * Variable is approximately feasible, i.e., every y[ o ][ l ] is binary
  * and sum_{ o , l } y[ o ][ l ] <= 1, both up to a tolerance eps. This
  * requires that generate_abstract_variables() has been called prior to
  * this method, otherwise false is returned. The tolerance eps is found as
  * follows:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< FNumber >, then
  *   it is fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration< FNumber >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0.
  *
  * The parameter useabstract is ignored, since the check is the same. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) optimal
 /** Returns true if the solution encoded in the current value of the y[][]
  * Variable is approximately feasible [see is_feasible()] and
  * approximately optimal for the current Objective. Since at most one y is
  * nonzero, the optimal value of the DiscreteSatelliteBlock is the constant
  * term of the Objective plus the minimum between 0 and the smallest
  * coefficient, which is compared with the current value of the Objective.
  * This requires that generate_abstract_variables() and
  * generate_objective() have been called prior to this method, otherwise
  * false is returned.
  *
  * The two tolerances, the relative one ceps for the objective value and
  * feps for the feasibility, are found as follows:
  *
  * - if optc is not nullptr and it is a
  *   SimpleConfiguration< std::pair< CNumber , FNumber > >, then
  *   ceps = optc->f_value.first and feps = optc->f_value.second;
  *
  * - if optc is not nullptr and it is a SimpleConfiguration< CNumber >, then
  *   ceps = optc->f_value, while feps is taken out of
  *   f_BlockConfig->f_is_feasible_Configuration as in is_feasible();
  *
  * - otherwise, if f_BlockConfig is not nullptr, then feps is taken
  *   out of f_BlockConfig->f_is_feasible_Configuration, while ceps
  *   is taken out of f_BlockConfig->f_is_optimal_Configuration
  *   assuming the latter is a SimpleConfiguration< CNumber >;
  *
  * - otherwise, ceps == feps == 0. */

 bool is_optimal( bool useabstract = false ,
                  Configuration * optc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a DiscreteSatelliteSolution with the current values of the y
 /** Returns a DiscreteSatelliteSolution sized for the OrbitSet x ell
  * Variable y[][] of this DiscreteSatelliteBlock; unless emptys == true, it
  * also holds their current values. The Configuration solc is ignored. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the DiscreteSatelliteBlock
 /** Method for handling Modification. Any "abstract Modification" (one for
  * which Modification::concerns_Block() is true) has its concerns_Block()
  * value set to false and is passed to guts_of_add_Modification() before
  * being forwarded to Block::add_Modification(). Since the
  * DiscreteSatelliteBlock has no "physical representation" besides the
  * size of its grid, the abstract Modification require no translation and
  * guts_of_add_Modification() does nothing; changes of the coefficients of
  * the Objective (e.g., those of a Lagrangian relaxation) are just
  * forwarded to the Solver. */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*-------- METHODS FOR PRINTING & SAVING THE DiscreteSatelliteBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the DiscreteSatelliteBlock
 *  @{ */

 /// print the DiscreteSatelliteBlock on an ostream
 /** Prints the size of the (orbit , level) grid of the
  * DiscreteSatelliteBlock, in the format read by load( std::istream & );
  * vlvl is ignored. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 unsigned char AR1; ///< nonzero if the Objective has been constructed
 unsigned char AR2; ///< nonzero if the Constraint have been constructed
 unsigned char AR3; ///< nonzero if the Variable have been constructed

 static constexpr unsigned char HasVar = 1;
 ///< bit of AR3 set if the Variable are constructed
 static constexpr unsigned char HasObj = 2;
 ///< bit of AR1 set if the Objective is constructed
 static constexpr unsigned char HasCnst = 4;
 ///< bit of AR2 set if the Constraint are constructed

 Index OrbitSet; ///< the number of candidate orbits
 Index ell;      ///< the number of levels of the threshold theta

 boost::multi_array< ColVariable , 2 > y; ///< the y[][] Variable

 FRowConstraint one; ///< the constraint sum_{ o , l } y[ o ][ l ] <= 1

 FRealObjective c; ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );

/*--------------------------------------------------------------------------*/
/*----------------------------- PRIVATE FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h; // insert DiscreteSatelliteBlock in the factory

/*--------------------------------------------------------------------------*/

 }; // end( class( DiscreteSatelliteBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS DiscreteSatelliteSolution ---------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a DiscreteSatelliteBlock
/** A DiscreteSatelliteSolution holds the values of the OrbitSet x ell
 * Variable y[][] of a DiscreteSatelliteBlock, in "row-major" order (that of
 * y[ o ][ l ] has index o * ell + l). An empty DiscreteSatelliteSolution
 * (one with no values) reads and writes nothing. */

class DiscreteSatelliteSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*-------------------------------- FRIENDS ---------------------------------*/

 friend DiscreteSatelliteBlock; ///< make DiscreteSatelliteBlock friend

/*--------- CONSTRUCTING AND DESTRUCTING DiscreteSatelliteSolution ---------*/

 /// constructor, it has nothing to do

 explicit DiscreteSatelliteSolution( void ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// de-serialize a DiscreteSatelliteSolution out of netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// destructor: it is virtual, and empty

 ~DiscreteSatelliteSolution() override = default;

/*----- METHODS DESCRIBING THE BEHAVIOR OF A DiscreteSatelliteSolution -----*/

 /// read the values of the y[][] out of the given DiscreteSatelliteBlock
 /** Reads the current values of the y[][] Variable of the given
  * DiscreteSatelliteBlock, which must have the same size as this
  * DiscreteSatelliteSolution (exception is thrown otherwise); does nothing
  * if the DiscreteSatelliteSolution is empty. */

 void read( const Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// write the values of the y[][] into the given DiscreteSatelliteBlock
 /** Writes the stored values into the y[][] Variable of the given
  * DiscreteSatelliteBlock, which must have the same size as this
  * DiscreteSatelliteSolution (exception is thrown otherwise); does nothing
  * if the DiscreteSatelliteSolution is empty. */

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// serialize a DiscreteSatelliteSolution into a netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a DiscreteSatelliteSolution with the values multiplied by factor

 DiscreteSatelliteSolution * scale( double factor ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// adds the given DiscreteSatelliteSolution, multiplied by multiplier
 /** Adds multiplier times the values of the given DiscreteSatelliteSolution,
  * which must have the same size, to the stored ones. */

 void sum( const Solution * solution , double multiplier ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a copy of the DiscreteSatelliteSolution
 /** Returns a copy of the DiscreteSatelliteSolution; if empty == true the
  * copy has the same size but all its values are zero. */

 DiscreteSatelliteSolution * clone( bool empty = false ) const override final;

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// print the DiscreteSatelliteSolution

 void print( std::ostream & output ) const override final;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*----------------------------- PRIVATE FIELDS -----------------------------*/

 DiscreteSatelliteBlock::Vec_FNumber v_y; ///< the values of the y[][]

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( DiscreteSatelliteSolution ) )

/** @} end( group( DiscreteSatelliteBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* DiscreteSatelliteBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File DiscreteSatelliteBlock.h --------------------*/
/*--------------------------------------------------------------------------*/
