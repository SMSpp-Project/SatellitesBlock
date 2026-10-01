/*--------------------------------------------------------------------------*/
/*----------------------- File MultiTargetBlockv2.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiTargetBlockv2, which implements
 * the Block concept [see Block.h] for the Satellite Constellation Design
 * Problem with the objective of the average maximum revisit time of the
 * targets, as a single Block, and for the class MultiTargetSolution
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

#ifndef __MultiTargetBlockv2
 #define __MultiTargetBlockv2
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
class MultiTargetBlockv2; // forward declaration of MultiTargetBlockv2

class MultiTargetSolution; // forward declaration of MultiTargetSolution

/*--------------------------------------------------------------------------*/
/*-------------------- MultiTargetBlockv2-RELATED TYPES --------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlockv2_TYPES MultiTargetBlockv2-related types
 *  @{ */

using p_MultiTargetBlockv2 = MultiTargetBlockv2 *;
///< a pointer to MultiTargetBlockv2

using Vec_MultiTargetBlockv2 = std::vector< p_MultiTargetBlockv2 >;
///< a vector of pointers to MultiTargetBlockv2

using Vec_MultiTargetBlockv2_it = Vec_MultiTargetBlockv2::iterator;
///< iterator for a Vec_MultiTargetBlockv2

using c_Vec_MultiTargetBlockv2 = const Vec_MultiTargetBlockv2;
///< a const vector of pointers to MultiTargetBlockv2

using c_Vec_MultiTargetBlockv2_it = c_Vec_MultiTargetBlockv2::iterator;
///< iterator for a c_Vec_MultiTargetBlockv2

/** @} end( group( MultiTargetBlockv2_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlockv2_CLASSES Classes in MultiTargetBlockv2.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MultiTargetBlockv2 ------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the multi-target constellation design problem, as a single Block
/** MultiTargetBlockv2 represents the same problem as MultiTargetBlock [see
 * MultiTargetBlock.h], i.e., minimizing the average maximum revisit time of
 * the targets with orbits common to all of them, but as a single Block
 * rather than as one SingleTargetBlock per target linked by "duplicate"
 * consistency constraints. The model is that of SingleTargetBlock [see
 * SingleTargetBlock.h], with the differences that
 *
 * - the orbit Variable activation[ i ][ c ] (satellite i, candidate orbit
 *   c) and the thresholds theta[ i ] are one single copy for all the
 *   targets, hence no "duplicate" constraints are needed;
 *
 * - every other Variable and Constraint of SingleTargetBlock has a further
 *   index, the last one, for the target, e.g., xi[ i ][ k ][ m ] is 1 if
 *   satellite i observes target m at time stamp k.
 *
 * MultiTargetBlockv2 has no sub-Block, and it is meant to be solved by a
 * general-purpose MILP Solver. In the comments of the fields, the numbers
 * of the constraints are those of SingleTargetBlock.h. */

class MultiTargetBlockv2 : public Block {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * MultiTargetBlockv2 defines three main public types:
 *
 * - FNumber, the type of the values of the Variable and of the data;
 *
 * - CNumber, the type of the objective coefficients;
 *
 * - FONumber, the type of the objective function value.
 *
 * They are all double; they are kept as separate names for uniformity with
 * the other Block of the module.
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

 friend MultiTargetSolution; ///< make MultiTargetSolution friend

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of MultiTargetBlockv2, taking a pointer to the father Block
 /** Constructor of MultiTargetBlockv2. It accepts a pointer to the father
  * Block, which can be of any type, defaulting to nullptr so that this can
  * also be used as the void constructor. */

 explicit MultiTargetBlockv2( Block * father = nullptr )
  : Block( father ) , AR1( 0 ) , AR2( 0 ) , AR3( 0 ) , n( 0 ) , t( 0 ) ,
    OrbitSet( 0 ) , targets( 0 ) , dt( 0 ) , T( 0 ) , alphaHalf( 0 ) ,
    thetaValFinal( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of MultiTargetBlockv2: deletes the abstract representation

 virtual ~MultiTargetBlockv2() { guts_of_destructor(); }

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
 /** Loads the instance out of an istream, which has the format of
  * MultiTargetBlock::load() [see there], and computes the candidate orbits
  * and the distances of their ground tracks from the targets in the same
  * way. Any previous abstract representation is deleted. Exception is
  * thrown if the input is not correct. If there is any Solver attached to
  * this MultiTargetBlockv2 then a NBModification (the "nuclear option") is
  * issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the MultiTargetBlockv2

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the MultiTargetBlockv2
 /** Generates the constraints of SingleTargetBlock, one copy per target
  * (save those involving only activation and theta, of which there is a
  * single copy). */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the Objective, the average maximum revisit time

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- Methods for reading the data of the MultiTargetBlockv2 ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the MultiTargetBlockv2
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of Deltat[ k ], the revisit time of target k

 FNumber get_Deltat( Index k ) const { return( Deltat[ k ].get_value() ); }

/*--------------------------------------------------------------------------*/
 /// sets the values of the Deltat in the range rng out of fstrt
 /** Sets the values of the Variable Deltat[ k ], for k in the range
  * [ rng.first , min( rng.second , number of targets ) ), to the
  * consecutive values starting from fstrt. */

 void set_Deltat( c_Vec_FNumber_it fstrt ,
                  Range rng = Range( 0 , Inf< Index >() ) );

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the
  * Variable of the MultiTargetBlockv2 is approximately feasible, i.e.,
  * every Variable is within its bounds and integer if it has to be, and
  * every static Constraint is satisfied, up to a tolerance eps; this
  * requires that the Variable and the Constraint have been generated,
  * otherwise false is returned. The tolerance eps is found as follows:
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
  * The tolerance is absolute for the Variable and relative for the
  * Constraint; the parameter useabstract is ignored. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// optimality of a MultiTargetBlockv2 is not checked: returns false
 /** Deciding whether the current solution is optimal would require to
  * solve the MultiTargetBlockv2, which is a MILP: this method always
  * returns false. */

 bool is_optimal( bool useabstract = false ,
                  Configuration * optc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a MultiTargetSolution with the current values of the Deltat
 /** Returns a MultiTargetSolution sized for the Variable Deltat[] of this
  * MultiTargetBlockv2; unless emptys == true, it also holds their current
  * values. The Configuration solc is ignored. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the MultiTargetBlockv2
 /** Method for handling Modification. Any "abstract Modification" (one for
  * which Modification::concerns_Block() is true) has its concerns_Block()
  * value set to false and is passed to guts_of_add_Modification() before
  * being forwarded to Block::add_Modification(). Since the data of the
  * MultiTargetBlockv2 are not changed by its abstract representation, the
  * abstract Modification require no translation and
  * guts_of_add_Modification() does nothing. */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE MultiTargetBlockv2 ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the MultiTargetBlockv2
 *  @{ */

 /// print the MultiTargetBlockv2 on an ostream
 /** Prints the size of the MultiTargetBlockv2 (targets, satellites, time
  * stamps, candidate orbits); vlvl is ignored. */

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

 Index n;        ///< the number of satellites
 Index t;        ///< the number of time stamps
 Index OrbitSet; ///< the number of candidate orbits
 Index targets;  ///< the number of targets

 FNumber dt;            ///< the time step [s]
 FNumber T;             ///< the time horizon [s]
 FNumber alphaHalf;     ///< the half-aperture of the sensors
 FNumber thetaValFinal; ///< the maximum threshold theta^{\max}

 boost::multi_array< double , 3 > CoverageSatLat;
 ///< CoverageSatLat[ m ][ k ][ c ], the latitude distances
 boost::multi_array< double , 3 > CoverageSatLong;
 ///< CoverageSatLong[ m ][ k ][ c ], the longitude distances

 std::vector< ColVariable > theta;  ///< the thresholds of the satellites
 std::vector< ColVariable > Deltat; ///< the maximum revisit times
 boost::multi_array< ColVariable , 2 > Deltat_k1;
 ///< the forward revisit times
 boost::multi_array< ColVariable , 2 > Deltat_k2;
 ///< the backward revisit times
 boost::multi_array< ColVariable , 2 > Deltat_k1A;
 ///< auxiliary, one per pair k < j
 boost::multi_array< ColVariable , 2 > Deltat_k2A;
 ///< auxiliary, one per pair k < j
 boost::multi_array< ColVariable , 2 > zeta; ///< the observations z
 boost::multi_array< ColVariable , 2 > b1;   ///< the forward indicators
 boost::multi_array< ColVariable , 2 > b2;   ///< the backward indicators
 boost::multi_array< ColVariable , 2 > d1;   ///< the forward min indicators
 boost::multi_array< ColVariable , 2 > d2;   ///< the backward min indicators
 boost::multi_array< ColVariable , 2 > h;    ///< the pairs h, k < j
 boost::multi_array< ColVariable , 3 > xi;
 ///< the observations xi[ i ][ k ][ m ]
 boost::multi_array< ColVariable , 2 > activation;
 ///< the orbit Variable pi[ i ][ c ]

 std::vector< FRowConstraint > orbitSelection; ///< sum_c pi[ i ][ c ] = 1
 std::vector< FRowConstraint > theta_UB;       ///< theta[ i ] <= theta^{\max}

 std::vector< FRowConstraint > Deltat_max_dt1; ///< Deltat >= dt
 boost::multi_array< FRowConstraint , 2 > Deltat_max1;
 ///< the forward part of (4)
 boost::multi_array< FRowConstraint , 2 > Deltat_max11;
 ///< b1 indicators of (4)
 boost::multi_array< FRowConstraint , 2 > Deltat_max2;
 ///< the backward part of (4)
 boost::multi_array< FRowConstraint , 2 > Deltat_max22;
 ///< b2 indicators of (4)

 boost::multi_array< FRowConstraint , 2 > Deltat_min_k1_1;
 ///< the linearized (2)
 boost::multi_array< FRowConstraint , 2 > Deltat_min_k1_2;
 ///< the linearized (2)
 boost::multi_array< FRowConstraint , 2 > d1_cnst;
 ///< one d1 per time stamp in (2)
 boost::multi_array< FRowConstraint , 2 > Deltat_min_k2_1;
 ///< the linearized (3)
 boost::multi_array< FRowConstraint , 2 > Deltat_min_k2_2;
 ///< the linearized (3)
 boost::multi_array< FRowConstraint , 2 > d2_cnst;
 ///< one d2 per time stamp in (3)

 boost::multi_array< FRowConstraint , 3 > activationSat_cnst;
 ///< xi <= z in (6)
 boost::multi_array< FRowConstraint , 2 > activationSat1_cnst;
 ///< z <= sum_i xi in (6)
 boost::multi_array< FRowConstraint , 2 > observation1;
 ///< sum_i xi <= 1 in (6)

 boost::multi_array< FRowConstraint , 2 > h_cnst_1; ///< the first of (7)
 boost::multi_array< FRowConstraint , 2 > h_cnst_2; ///< the second of (7)
 boost::multi_array< FRowConstraint , 2 > h_cnst_3; ///< the third of (7)

 boost::multi_array< FRowConstraint , 3 > obs2_cnst;
 ///< the constraints (5) for the latitude
 boost::multi_array< FRowConstraint , 3 > obs4_cnst;
 ///< the constraints (5) for the longitude

 std::vector< FRowConstraint > obs_cnst_h; ///< sum h >= 1, per target
 boost::multi_array< FRowConstraint , 2 > obs_cnst_xi;
 ///< sum_k xi >= 3, per target

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

 SMSpp_insert_in_factory_h; // insert MultiTargetBlockv2 in the factory

/*--------------------------------------------------------------------------*/

 }; // end( class( MultiTargetBlockv2 ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS MultiTargetSolution ------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a MultiTargetBlockv2
/** A MultiTargetSolution holds the values of the Variable Deltat[] of a
 * MultiTargetBlockv2, i.e., the maximum revisit times of the targets. An
 * empty MultiTargetSolution (one with no values) reads and writes
 * nothing. */

class MultiTargetSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*-------------------------------- FRIENDS ---------------------------------*/

 friend MultiTargetBlockv2; ///< make MultiTargetBlockv2 friend

/*------------ CONSTRUCTING AND DESTRUCTING MultiTargetSolution ------------*/

 /// constructor, it has nothing to do

 explicit MultiTargetSolution( void ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// de-serialize a MultiTargetSolution out of netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// destructor: it is virtual, and empty

 ~MultiTargetSolution() override = default;

/*-------- METHODS DESCRIBING THE BEHAVIOR OF A MultiTargetSolution --------*/

 /// read the values of the Deltat out of the given MultiTargetBlockv2
 /** Reads the current values of the Variable Deltat[] of the given
  * MultiTargetBlockv2, which must have as many targets as the values of
  * this MultiTargetSolution (exception is thrown otherwise); does nothing
  * if the MultiTargetSolution is empty. */

 void read( const Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// write the values of the Deltat into the given MultiTargetBlockv2
 /** Writes the stored values into the Variable Deltat[] of the given
  * MultiTargetBlockv2; does nothing if the MultiTargetSolution is empty. */

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// serialize a MultiTargetSolution into a netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a MultiTargetSolution with the values multiplied by factor

 MultiTargetSolution * scale( double factor ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// adds the given MultiTargetSolution, multiplied by multiplier
 /** Adds multiplier times the values of the given MultiTargetSolution,
  * which must have the same size, to the stored ones. */

 void sum( const Solution * solution , double multiplier ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a copy of the MultiTargetSolution
 /** Returns a copy of the MultiTargetSolution; if empty == true the copy
  * has the same size but all its values are zero. */

 MultiTargetSolution * clone( bool empty = false ) const override final;

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// print the MultiTargetSolution

 void print( std::ostream & output ) const override final;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*----------------------------- PRIVATE FIELDS -----------------------------*/

 MultiTargetBlockv2::Vec_FNumber v_Deltat; ///< the values of the Deltat

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( MultiTargetSolution ) )

/** @} end( group( MultiTargetBlockv2_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* MultiTargetBlockv2.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File MultiTargetBlockv2.h ----------------------*/
/*--------------------------------------------------------------------------*/
