/*--------------------------------------------------------------------------*/
/*------------------------ File SingleTargetBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SingleTargetBlock, which implements
 * the Block concept [see Block.h] for the revisit-time problem of a single
 * target of the Satellite Constellation Design Problem, and for the class
 * SingleTargetSolution holding a solution of it.
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

#ifndef __SingleTargetBlock
 #define __SingleTargetBlock
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
class SingleTargetBlock; // forward declaration of SingleTargetBlock

class SingleTargetSolution; // forward declaration of SingleTargetSolution

/*--------------------------------------------------------------------------*/
/*-------------------- SingleTargetBlock-RELATED TYPES ---------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleTargetBlock_TYPES SingleTargetBlock-related types
 *  @{ */

using p_SingleTargetBlock = SingleTargetBlock *;
///< a pointer to SingleTargetBlock

using Vec_SingleTargetBlock = std::vector< p_SingleTargetBlock >;
///< a vector of pointers to SingleTargetBlock

using Vec_SingleTargetBlock_it = Vec_SingleTargetBlock::iterator;
///< iterator for a Vec_SingleTargetBlock

using c_Vec_SingleTargetBlock = const Vec_SingleTargetBlock;
///< a const vector of pointers to SingleTargetBlock

using c_Vec_SingleTargetBlock_it = c_Vec_SingleTargetBlock::iterator;
///< iterator for a c_Vec_SingleTargetBlock

/** @} end( group( SingleTargetBlock_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleTargetBlock_CLASSES Classes in SingleTargetBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SingleTargetBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the revisit-time problem of a single target
/** SingleTargetBlock represents the problem of choosing the orbits of a
 * constellation of s satellites so as to minimize the maximum revisit time
 * of a single target, i.e., the longest time between two consecutive
 * observations of the target; it is the sub-Block, one per target, of a
 * MultiTargetBlock [see MultiTargetBlock.h]. Let T be the time horizon, dt
 * the time step, [t] = { 0 , ... , t - 1 } the set of the t = T / dt time
 * stamps and [C] the set of the candidate orbits, common to all the
 * satellites, for each of which the latitude and longitude distances
 * \f$ \Delta lat[ k ][ c ] \f$ and \f$ \Delta long[ k ][ c ] \f$ between
 * the ground track and the target at time stamp k are known (they are
 * computed by MultiTargetBlock::load()). The binary Variable are
 *
 * - \f$ \pi[ i ][ c ] \f$ (activation), 1 if satellite i is placed on orbit
 *   c;
 *
 * - \f$ \xi[ i ][ k ] \f$ (xi), 1 if satellite i observes the target at
 *   time stamp k;
 *
 * - \f$ z[ k ] \f$ (zeta), 1 if some satellite observes the target at time
 *   stamp k;
 *
 * - \f$ h[ k , j ] \f$ (h), for k < j, 1 if the target is observed at both
 *   time stamps k and j;
 *
 * - b1, b2, d1, d2, the indicators used to linearize the min and max
 *   operators below;
 *
 * and the continuous ones are the thresholds \f$ \theta[ i ] \f$ of the
 * satellites, the maximum revisit time \f$ \Delta \f$ (Deltat) and the
 * forward and backward revisit times \f$ \Delta^1[ k ] \f$ and
 * \f$ \Delta^2[ k ] \f$ (Deltat_k1 and Deltat_k2) at every time stamp k.
 * The model is
 *
 * \f[
 *  \min \; \Delta / m                                             \qquad (1)
 * \f]
 * \f[
 *  \Delta^1[ k ] = \min_{ j > k } \{ dt ( j - k ) h[ k , j ] +
 *  T / 2 ( 1 - h[ k , j ] ) \}                                    \qquad (2)
 * \f]
 * \f[
 *  \Delta^2[ k ] = \min_{ j > k } \{ ( T - dt ( j - k ) ) h[ k , j ] +
 *  T / 2 ( 1 - h[ k , j ] ) \}                                    \qquad (3)
 * \f]
 * \f[
 *  \Delta \geq dt \; , \quad \Delta \geq \Delta^1[ k ] - T / 2 \, b1[ k ]
 *  \; , \quad \Delta \geq \Delta^2[ k ] - T / 2 \, b2[ k ]        \qquad (4)
 * \f]
 * \f[
 *  \xi[ i ][ k ] = 1 \Longrightarrow \max \Bigl\{ \sum_{ c \in [C] }
 *  \pi[ i ][ c ] \Delta lat[ k ][ c ] \, , \, \sum_{ c \in [C] }
 *  \pi[ i ][ c ] \Delta long[ k ][ c ] \Bigr\} \leq \theta[ i ]   \qquad (5)
 * \f]
 * \f[
 *  \xi[ i ][ k ] \leq z[ k ] \leq \sum_{ i } \xi[ i ][ k ] \; , \quad
 *  \sum_{ i } \xi[ i ][ k ] \leq 1                                \qquad (6)
 * \f]
 * \f[
 *  h[ k , j ] \leq z[ k ] \; , \quad h[ k , j ] \leq z[ j ] \; , \quad
 *  h[ k , j ] \geq z[ k ] + z[ j ] - 1                            \qquad (7)
 * \f]
 *
 * together with \f$ \sum_{ c } \pi[ i ][ c ] = 1 \f$, \f$ \theta[ i ] \leq
 * \theta^{\max} \f$, \f$ \sum_{ k < j } h[ k , j ] \geq 1 \f$ (the target is
 * observed at least twice) and \f$ \sum_k \xi[ i ][ k ] \geq 3 \f$ (every
 * satellite observes the target at least three times). The objective (1)
 * is the maximum revisit time divided by the number m of targets, so that
 * the objective of a MultiTargetBlock, the sum of those of its
 * SingleTargetBlock, is the average maximum revisit time. Constraints (2)
 * and (3) are linearized with the indicators d1 and d2, (4) with b1 and
 * b2, and (5) with big-M constraints, M being the largest distance over
 * [C]; (7) is the linearization of \f$ h[ k , j ] = z[ k ] z[ j ] \f$. */

class SingleTargetBlock : public Block {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * SingleTargetBlock defines three main public types:
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

 friend SingleTargetSolution; ///< make SingleTargetSolution friend

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of SingleTargetBlock, taking a pointer to the father Block
 /** Constructor of SingleTargetBlock. It accepts a pointer to the father
  * Block, which can be of any type, defaulting to nullptr so that this can
  * also be used as the void constructor. */

 explicit SingleTargetBlock( Block * father = nullptr )
  : Block( father ) , AR1( 0 ) , AR2( 0 ) , AR3( 0 ) , n( 0 ) , t( 0 ) ,
    OrbitSet( 0 ) , targets( 0 ) , dt( 0 ) , T( 0 ) , alphaHalf( 0 ) ,
    thetaVal( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of SingleTargetBlock: deletes the abstract representation

 virtual ~SingleTargetBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance out of memory
 /** Loads the instance out of memory; the data is computed by
  * MultiTargetBlock::load(), which reads them out of a file. The
  * parameters are:
  *
  * - num_targets, the number m of targets of the MultiTargetBlock (the
  *   objective is divided by it);
  *
  * - num_satellites, the number s of satellites;
  *
  * - time_step, the time step dt [s];
  *
  * - horizon, the time horizon T [s], so that there are T / dt time
  *   stamps;
  *
  * - thetaValues, the maximum threshold \f$ \theta^{\max} \f$;
  *
  * - indOrbit, the number of candidate orbits [C];
  *
  * - aHalf, the half-aperture of the sensors;
  *
  * - CoverageLat and CoverageLong, with CoverageLat[ k ][ c ] the latitude
  *   distance between the ground track of orbit c and the target at time
  *   stamp k, and the same for the longitude.
  *
  * Any previous abstract representation is deleted. If there is any Solver
  * attached to this SingleTargetBlock then a NBModification (the "nuclear
  * option") is issued. */

 void load( Index num_targets , Index num_satellites , FNumber time_step ,
            FNumber horizon , FNumber thetaValues , Index indOrbit ,
            FNumber aHalf ,
            const boost::multi_array< double , 2 > & CoverageLat ,
            const boost::multi_array< double , 2 > & CoverageLong );

/*--------------------------------------------------------------------------*/
 /// loading a SingleTargetBlock out of an istream is not supported
 /** A SingleTargetBlock is only loaded by its MultiTargetBlock, out of
  * memory: this always throws std::logic_error. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// loads the instance from the file with the given name
 /** The version of Block, which opens the file and calls
  * load( std::istream & ). */

 using Block::load;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the SingleTargetBlock

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the SingleTargetBlock
 /** Generates the linearized constraints (2) - (7) of the class comments,
  * and the other ones listed there. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the Objective (1) of the SingleTargetBlock

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------- Methods for reading the data of the SingleTargetBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the SingleTargetBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the Variable activation[ i ][ c ]
 /** Returns a pointer to the Variable activation[ i ][ c ], i.e.,
  * \f$ \pi[ i ][ c ] \f$; this requires that generate_abstract_variables()
  * has been called. */

 [[nodiscard]] ColVariable * i2p_pi( Index i , Index c ) const {
  return( const_cast< ColVariable * >( &activation[ i ][ c ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the Variable theta[ i ]
 /** Returns a pointer to the threshold Variable theta[ i ] of satellite i;
  * this requires that generate_abstract_variables() has been called. */

 [[nodiscard]] ColVariable * i2p_theta( Index i ) const {
  return( const_cast< ColVariable * >( &theta[ i ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable activation[ i ][ c ]

 FNumber get_activation( Index i , Index c ) const {
  return( activation[ i ][ c ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable activation[ i ][ c ] and fixes it

 void set_activation( Index i , Index c , int value ) {
  activation[ i ][ c ].set_value( value );
  activation[ i ][ c ].is_fixed( true );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the maximum revisit time Deltat

 FNumber get_Deltat( void ) const { return( Deltat[ 0 ].get_value() ); }

/*--------------------------------------------------------------------------*/
 /// sets the value of Deltat out of the range rng of the values in fstrt
 /** Sets the value of the Variable Deltat to *fstrt if rng.first == 0
  * (Deltat being a single Variable, there is nothing to do otherwise). */

 void set_Deltat( c_Vec_FNumber_it fstrt ,
                  Range rng = Range( 0 , Inf< Index >() ) );

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the
  * Variable of the SingleTargetBlock is approximately feasible, i.e., every
  * Variable is within its bounds and integer if it has to be, and every
  * static Constraint is satisfied, up to a tolerance eps; this requires
  * that the Variable and the Constraint have been generated, otherwise
  * false is returned. The tolerance eps is found as follows:
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
 /// optimality of a SingleTargetBlock is not checked: returns false
 /** Deciding whether the current solution is optimal would require to
  * solve the SingleTargetBlock, which is a MILP: this method always
  * returns false. */

 bool is_optimal( bool useabstract = false ,
                  Configuration * optc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a SingleTargetSolution with the current value of Deltat
 /** Returns a SingleTargetSolution sized for the Variable Deltat of this
  * SingleTargetBlock; unless emptys == true, it also holds its current
  * value. The Configuration solc is ignored. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the SingleTargetBlock
 /** Method for handling Modification. Any "abstract Modification" (one for
  * which Modification::concerns_Block() is true) has its concerns_Block()
  * value set to false and is passed to guts_of_add_Modification() before
  * being forwarded to Block::add_Modification(). Since the data of the
  * SingleTargetBlock are not changed by its abstract representation, the
  * abstract Modification require no translation and
  * guts_of_add_Modification() does nothing. */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE SingleTargetBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SingleTargetBlock
 *  @{ */

 /// print the SingleTargetBlock on an ostream
 /** Prints the size of the SingleTargetBlock (satellites, time stamps,
  * candidate orbits) and its maximum threshold; vlvl is ignored. */

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
 Index targets;  ///< the number of targets of the MultiTargetBlock

 FNumber dt;        ///< the time step [s]
 FNumber T;         ///< the time horizon [s]
 FNumber alphaHalf; ///< the half-aperture of the sensors
 FNumber thetaVal;  ///< the maximum threshold theta^{\max}

 boost::multi_array< double , 2 > CoverageSatLat;
 ///< CoverageSatLat[ k ][ c ], the latitude distances
 boost::multi_array< double , 2 > CoverageSatLong;
 ///< CoverageSatLong[ k ][ c ], the longitude distances

 std::vector< ColVariable > Deltat;     ///< the maximum revisit time
 std::vector< ColVariable > Deltat_k1;  ///< the forward revisit times (2)
 std::vector< ColVariable > Deltat_k2;  ///< the backward revisit times (3)
 std::vector< ColVariable > Deltat_k1A; ///< auxiliary, one per pair k < j
 std::vector< ColVariable > Deltat_k2A; ///< auxiliary, one per pair k < j
 std::vector< ColVariable > zeta;       ///< the observations z[ k ]
 std::vector< ColVariable > b1;         ///< the indicators of (4), forward
 std::vector< ColVariable > b2;         ///< the indicators of (4), backward
 std::vector< ColVariable > d1;         ///< the indicators of the min (2)
 std::vector< ColVariable > d2;         ///< the indicators of the min (3)
 std::vector< ColVariable > h;          ///< the pairs h[ k , j ], k < j
 std::vector< ColVariable > theta;      ///< the thresholds of the satellites

 boost::multi_array< ColVariable , 2 > xi;
 ///< the observation Variable xi[ i ][ k ]
 boost::multi_array< ColVariable , 2 > activation;
 ///< the orbit Variable pi[ i ][ c ]

 std::vector< FRowConstraint > orbitSelection; ///< sum_c pi[ i ][ c ] = 1
 std::vector< FRowConstraint > theta_UB;       ///< theta[ i ] <= theta^{\max}

 std::vector< FRowConstraint > Deltat_max_dt1; ///< Deltat >= dt
 std::vector< FRowConstraint > Deltat_max1;    ///< the forward part of (4)
 std::vector< FRowConstraint > Deltat_max11;   ///< b1 indicators of (4)
 std::vector< FRowConstraint > Deltat_max2;    ///< the backward part of (4)
 std::vector< FRowConstraint > Deltat_max22;   ///< b2 indicators of (4)

 std::vector< FRowConstraint > Deltat_min_k1_1; ///< the linearized (2)
 std::vector< FRowConstraint > Deltat_min_k1_2; ///< the linearized (2)
 std::vector< FRowConstraint > d1_cnst; ///< one d1 per time stamp in (2)
 std::vector< FRowConstraint > Deltat_min_k2_1; ///< the linearized (3)
 std::vector< FRowConstraint > Deltat_min_k2_2; ///< the linearized (3)
 std::vector< FRowConstraint > d2_cnst; ///< one d2 per time stamp in (3)

 boost::multi_array< FRowConstraint , 2 > activationSat_cnst;
 ///< xi[ i ][ k ] <= z[ k ] in (6)
 std::vector< FRowConstraint > activationSat1_cnst;
 ///< z[ k ] <= sum_i xi[ i ][ k ] in (6)
 std::vector< FRowConstraint > observation1;
 ///< sum_i xi[ i ][ k ] <= 1 in (6)

 std::vector< FRowConstraint > h_cnst_1; ///< h[ k , j ] <= z[ k ] in (7)
 std::vector< FRowConstraint > h_cnst_2; ///< h[ k , j ] <= z[ j ] in (7)
 std::vector< FRowConstraint > h_cnst_3;
 ///< h[ k , j ] >= z[ k ] + z[ j ] - 1 in (7)

 boost::multi_array< FRowConstraint , 2 > obs2_cnst;
 ///< the constraints (5) for the latitude
 boost::multi_array< FRowConstraint , 2 > obs4_cnst;
 ///< the constraints (5) for the longitude

 std::vector< FRowConstraint > obs_cnst_h;  ///< sum h >= 1
 std::vector< FRowConstraint > obs_cnst_xi; ///< sum_k xi[ i ][ k ] >= 3

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

 SMSpp_insert_in_factory_h; // insert SingleTargetBlock in the Block factory

/*--------------------------------------------------------------------------*/

 }; // end( class( SingleTargetBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SingleTargetSolution -----------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a SingleTargetBlock
/** A SingleTargetSolution holds the value of the Variable Deltat of a
 * SingleTargetBlock, i.e., the maximum revisit time of the target. An
 * empty SingleTargetSolution (one with no value) reads and writes
 * nothing. */

class SingleTargetSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*-------------------------------- FRIENDS ---------------------------------*/

 friend SingleTargetBlock; ///< make SingleTargetBlock friend

/*----------- CONSTRUCTING AND DESTRUCTING SingleTargetSolution ------------*/

 /// constructor, it has nothing to do

 explicit SingleTargetSolution( void ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// de-serialize a SingleTargetSolution out of netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// destructor: it is virtual, and empty

 ~SingleTargetSolution() override = default;

/*------- METHODS DESCRIBING THE BEHAVIOR OF A SingleTargetSolution --------*/

 /// read the value of Deltat out of the given SingleTargetBlock
 /** Reads the current value of the Variable Deltat of the given
  * SingleTargetBlock; does nothing if the SingleTargetSolution is empty. */

 void read( const Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// write the value of Deltat into the given SingleTargetBlock
 /** Writes the stored value into the Variable Deltat of the given
  * SingleTargetBlock; does nothing if the SingleTargetSolution is empty. */

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// serialize a SingleTargetSolution into a netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a SingleTargetSolution with the value multiplied by factor

 SingleTargetSolution * scale( double factor ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// adds the given SingleTargetSolution, multiplied by multiplier
 /** Adds multiplier times the value of the given SingleTargetSolution,
  * which must have the same size, to the stored one. */

 void sum( const Solution * solution , double multiplier ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a copy of the SingleTargetSolution
 /** Returns a copy of the SingleTargetSolution; if empty == true the copy
  * has the same size but its value is zero. */

 SingleTargetSolution * clone( bool empty = false ) const override final;

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// print the SingleTargetSolution

 void print( std::ostream & output ) const override final;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*----------------------------- PRIVATE FIELDS -----------------------------*/

 SingleTargetBlock::Vec_FNumber v_Deltat; ///< the value of Deltat, if any

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( SingleTargetSolution ) )

/** @} end( group( SingleTargetBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SingleTargetBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File SingleTargetBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
