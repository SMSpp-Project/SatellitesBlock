/*--------------------------------------------------------------------------*/
/*-------------------- File SatelliteBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SatelliteBlock, which implements
 * the Block concept [see Block.h] for the solution of Satellite observavility 
 * problem.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SatelliteBlock
 #define __SatelliteBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class SatelliteBlock;     // forward declaration of SatelliteBlock

 class SatelliteSolution;  // forward declaration of SatelliteSolution

/*--------------------------------------------------------------------------*/
/*----------------------- SatelliteBlock-RELATED TYPES ---------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteBlock_TYPES SatelliteBlock-related types
 *  @{ */

 using p_SatelliteBlock = SatelliteBlock *;  ///< a pointer to SatelliteBlock

 using Vec_SatelliteBlock = std::vector< p_SatelliteBlock>;
 ///< a vector of pointers to SatelliteBlock

 using Vec_SatelliteBlock_it = Vec_SatelliteBlock::iterator;
 ///< iterator for a Vec_SatelliteBlock

 using c_Vec_SatelliteBlock = const Vec_SatelliteBlock;
 ///< a const vector of pointers to SatelliteBlock

 using c_Vec_SatelliteBlock_it = c_Vec_SatelliteBlock::iterator;
 ///< iterator for a c_Vec_SatelliteBlock

/** @}  end( group( SatelliteBlock_TYPES ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteBlock_CLASSES Classes in SatelliteBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SatelliteBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Implementation of a simple SatelliteBlock concept.
/*
* Let T and \Delta t_m be the time horizon and the time width between two 
* consecutive target observations for the target m\in\mathcal X$, respectively. 
* We indicate with [T] := \{ 1,2,...,T \} the set of time stamps. 
* [C] is the set of possible orbital configurations.
* We consider three sets of binary decision variables:
*
* (1) \pi[ c ] \in \{ 0 , 1 \}, c \in [C], indicating which orbit configuration c 
* is selected for the satellite, 
* (2) \xi[ t ][ m ] \in \{ 0 , 1 \} indicating whether the satellite observes the 
* target m at time t, and
* (3) zeta \in \{ 0 , 1 \} indicating whether the current satellite is active in 
* the constellation
* 
* Mathematically speaking, we have that the following constraints hold for each 
* satellite.
*
* \f[
* \xi[ t ][ m ]=0 \Longrightarrow \max \{ \sum_{c\in [C]} \pi[ c ]\, 
\Delta lat[ c ][ t ][ m ]}, \sum_{c\in [C]} \pi[ c ] \,\Delta long[ c ][ t ][ m ] \} 
\ge \theta^{\max}, \forall t\in T(dt), \forall m \in \mathcal{X}           (1)
* \f]
* \f[
* \sum_{ c \in [C] } \pi[ c ] = 1                                          (2)
* \f]
* \f[
* \xi[ t ][ m ] \leq zeta, \forall t\in T(dt), \forall m \in \mathcal{X}   (3)
* \f]
*
* Constraints (1) impose that, if the selected distances in longitude and latitude, 
* \Delta lat[ c ][ t ][ m ] and \Delta long[ c ][ t ][ m ], are smaller than the 
* threshold $\theta^{\max}$, then the target m is observed by the current satellite 
* at time t. Constraint (1) is opportunely linearized when defining the constraints
* of the current Block. Constraint (2) requires that exactly one configuration is 
* selected for the current satellite. Finally, constraints (3) active the current 
* satellite in the constellation if it observes at least one target.
*/

class SatelliteBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * SatelliteBlock defines three main public types:
 *
 * - FNumber, the type of flow variables, arc capacities, and node deficits;
 *
 * - CNumber, the type of flow costs, node potentials, and arc reduced costs;
 *
 * - FONumber, the type of objective function value.
 *
 * By re-defining the types in this section, some (but not all) solution
 * algorithms may be able to work with the "smallest" choice of data type 
 * that is capable of properly representing the data of the instances to be
 * solved. This may be relevant due to an important property of DCR problems:
 * *if all arc capacities and node deficits are integer, then there exists an
 * integral optimal primal solution*, and *if all arc costs are integer,
 * then there exists an integral optimal dual solution*. Even more
 * importantly, *many solution algorithms will in fact produce an integral
 * primal/dual solution for free*, because *every primal/dual solution they
 * generate during the solution process is naturally integral*. Therefore,
 * one can use integer data types to represent everything connected with
 * flows and/or costs if the corresponding data is integer in all instances
 * one needs to solve. This directly translates in significant memory savings
 * and/or speed improvements.
 *
 * However, while using a SatelliteBlock as a part of some larger problem, it may
 * be difficult to fully exploit this property: even if some Solver can
 * exploit it, not all of them may be able to (one example are Interior-Point
 * approaches, which require both flow and cost variables to be continuous),
 * and maybe some other aspects of the overall solution algorithm will require
 * general double data anyway. One should actually have Block template over
 * all these types to be able to fully exploit this property, which may be a
 * future evolution but is not what this implementation does. The current
 * choice is to use the "worst case scenario" where FNumber == CNumber ==
 * OFNumber == double, although the data types are left there and it is
 * therefore in principle possible to change this. Note, however, that the
 * above integrality property only holds for *linear* DCR problems. Should
 * the class be extended, by even allowing arc costs to be convex quadratic
 * (the simplest possible nonlinear extension), then a single arc with a
 * nonzero quadratic cost coefficient implies that optimal flows and
 * potentials may be fractional even if all the data of the problem
 * (comprised quadratic cost coefficients) is integer. Hence, for such a
 * setting FNumber == CNumber == OFNumber == double is actually *mandatory*,
 * for any reasonable algorithm will typically misbehave otherwise.
 @{ */

/*--------------------------------------------------------------------------*/

 typedef double FNumber;                     ///< type of arc flow / deficit
 typedef const FNumber c_FNumber;            ///< a read-only FNumber

 typedef std::vector< FNumber > Vec_FNumber; ///< a vector of FNumber
 typedef const Vec_FNumber c_Vec_FNumber;    ///< a const vector of FNumber

 typedef Vec_FNumber::iterator Vec_FNumber_it;   ///< iterator in Vec_FNumber
 typedef Vec_FNumber::const_iterator c_Vec_FNumber_it;
                                           ///< const iterator in Vec_FNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double CNumber;                     ///< type of arc cost / potential
 typedef const CNumber c_CNumber;            ///< a read-only CNumber

 typedef std::vector< CNumber > Vec_CNumber;  ///< a vector of CNumber
 typedef const Vec_CNumber c_Vec_CNumber;     ///< a const vector of CNumber

 typedef Vec_CNumber::iterator Vec_CNumber_it;   ///< iterator in Vec_CNumber
 typedef Vec_CNumber::const_iterator c_Vec_CNumber_it;
                                           ///< const iterator in Vec_CNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double FONumber; 
 /**< type of the objective function: has to hold sums of products of
    FNumber(s) by CNumber(s) */

 typedef const FONumber c_FONumber;             ///< a read-only FONumber

 typedef std::vector< FONumber > Vec_FONumber;  ///< a vector of FONumber
 typedef const Vec_FONumber c_Vec_FONumber;     ///< a const vector of FONumber

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend SatelliteSolution;  ///< make SatelliteSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of SatelliteBlock, taking a pointer to the father (generic) Block
 /** Constructor of SatelliteBlock. It accepts a pointer to the father Block, which
  * can be of any type, defaulting to nullptr so that this can also be used as
  * the void constructor. */

 explicit SatelliteBlock( Block *father = nullptr )
  : Block( father ) , AR1(0), AR2(0), AR3(0) { }
               

/*--------------------------------------------------------------------------*/
 /// destructor of SatelliteBlock: deletes the abstract representation, if any

 virtual ~SatelliteBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
  /*
  * Like load( std::istream & ), if there is any Solver attached to this
  * SatelliteBlock then a NBModification (the "nuclear option") is issued. */

 void load( FNumber n , FNumber dt , FNumber T , FNumber altValues , FNumber thetaValues ,
              FNumber indexOrbit, FNumber aHalf ,
              boost::multi_array< double , 3 > CoverageLat , boost::multi_array< double , 3 > CoverageLong );

 void load( std::istream &input , char frmt = 0 ) override;

 [[nodiscard]] ColVariable * i2p_r( Index iii , Index jjj ) const {
  return( const_cast< ColVariable * >( &xi[ iii ][ jjj ] ) );
 }

 [[nodiscard]] ColVariable * i2p_z() const {
  return( const_cast< ColVariable * >( &zeta[ 0 ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the SatelliteBlock
 /** Method that generates the abstract Variable of the Satellite. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the Satellite
 /** Method that generates the abstract constraint of the Satellite. */
 
 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the Satellite
 /** Method that generates the objective of the Satellite. */

 void generate_objective( Configuration *objc = nullptr ) override;

 //void generate_dynamic_constraints( Configuration *stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SatelliteBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the SatelliteBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }
  
/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SatelliteBlock is approximately feasible. This clearly
  * requires the Variable of the SatelliteBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method.
  *
  * The parameter for deciding what "approximately feasible" exactly means is
  * a single FNumber value, representing the *relative* tolerance for
  * satisfaction of both flow conservation constraint and flow upper/lower
  * bounds. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< FNumber >, then
  *   it is fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration< FNumber >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0. */
 
 bool is_feasible( bool useabstract = false , Configuration *fsbc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) optimal
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the SatelliteBlock is approximately optimal, which means that
  * it is approximately feasible, that the dual solution encoded in the
  * current value of the dual multipliers of both the flow conservation and
  * bound constraints is approximately feasible, and that the two
  * approximately satisfies the Complementary Slackness Conditions. This
  * clearly requires that both the Variable and the Constraint of the
  * SatelliteBlock to have been defined, i.e., that generate_abstract_variables()
  * and generate_abstract_constraints() have been called prior to this method.
  *
  * This requires two parameters for deciding what "approximately feasible"
  * means, one for the primal (feps) and one for the dual (ceps), like in
  * complementary_slackness(). These are found as follows:
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
 
 bool is_optimal( bool useabstract = false  , Configuration *optc = nullptr )
  override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the zeta variable

 FNumber get_zeta() const {
   return( zeta[ 0 ].get_value() );
  }

  void set_zeta( c_Vec_FNumber_it fstrt ,
		 Range rng = Range( 0 , Inf< Index >() ) );
  
 /// returns a SatelliteSolution representing the current solution of this SatelliteBlock

 Solution * get_Solution( Configuration *solc = nullptr ,
 			  bool emptys = true ) override;


/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// returns true if there is any Solver "listening to this SatelliteBlock"
 /** Returns true if there is any Solver "listening to this SatelliteBlock", or if
  * the SatelliteBlock has to "listen" anyway because the "abstract" representation
  * is constructed, and therefore "abstract" Modification have to be generated
  * anyway to keep the two representations in sync.
  *
  * No, this should not be needed. In fact, if the "abstract" representation
  * is modified with the default eModBlck value of issueMod, it is issued
  * irrespectively to the value of anyone_there(); see Observer::issue_mod().
  * If the value of issueMod is anything else the  "abstract" representation
  * has been modified already and there is no point in issuing the
  * Modification.
  * Note that that Observer::issue_mod() does not check if the "abstract"
  * representation has been constructed, but this is clearly not
  * necessary, as the Modification we are speaking of are issued while
  * changing the "abstract" representation, if that has not been
  * constructed then it cannot issue Modification

 bool anyone_there( void ) const override {
  return( AR ? true : Block::anyone_there() );
  }
 */
/*--------------------------------------------------------------------------*/
 /// adding a new Modification to the SatelliteBlock
 /** Method for handling Modification.
  *
  * The version of SatelliteBlock has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the SatelliteBlock, and "translate"
  * them into both changes of the actual data structures and corresponding
  * "physical Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the SatelliteBlock, and this must not be
  * repeated. In particular, this would be an issue if the Modification would
  * be [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * SatelliteBlock), which would mean that the "physical Modification" would be
  * issued twice.
  *
  * The following "abstract Modification" are handled:
  *
  * - GroupModification, that are simply unpacked into the individual
  *   sub-[Group]Modification and dealt with individually;
  *
  * - C05FunctionModRngd and C05FunctionModSbst changing coefficients coming
  *   from the (LinearFunction into the FRow)Objective, but *not* from the
  *   (LinearFunction into the FRow)Constraint;
  *
  * - RowConstraintMod changing the RHS of the bound constraints and both
  *   sides at once of the flow conservation ones, but not any other
  *   combination; and note that the RHS of the bound constraints may not
  *   be changeable at all if they have not been constructed, in which
  *   case there cannot be any Modification to handle here;
  *
  * - VariableMod fixing and un-fixing a flow ColVariable; however, note
  *   that *fixing is only permitted if the value() of the ColVariable is
  *   zero*, because that corresponds to closing the arc, exception being
  *   thrown otherwise.
  *
  * Any other Modification reaching the SatelliteBlock will lead to exception
  * being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE SatelliteBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SatelliteBlock
 *  @{ */

 /// print the SatelliteBlock on an ostream with the given verbosity
 /** Protected method to print information about the SatelliteBlock; with the
  * "complete" level ('C') it outputs the SatelliteBlock in DIMACS format. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 unsigned char AR1,AR2,AR3;               ///< bit-wise coded: what abstract is there
 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variable have been constructed
 static constexpr unsigned char HasObj = 2;
 ///< first bit of AR == 1 if the Objective have been constructed
 static constexpr unsigned char HasCnst = 3;
 ///< first bit of AR == 1 if the Constraint have been constructed

 Index n;                        ///< the number of targets
 Index t;                        ///< the total number of time step
 FNumber OrbitSet;               ///< the number of configurations

 FNumber indexOrbit;
 FNumber dt;                     ///< the time discretization step
 FNumber T;                      ///< the simulation time horizon
 FNumber alphaHalf;              ///< the alpha half satellite parameter

 FNumber altitudeVal;            ///< vector of satellite altitude values
 FNumber thetaVal;               ///< vector of satellite theta values

 boost::multi_array< double , 3 > CoverageSatLat;  ///< the matrix of pre-computed CoverageSatLat
 boost::multi_array< double , 3 > CoverageSatLong; ///< the matrix of pre-computed CoverageSatLong

 double f_cond_lower;            ///< conditional lower bound, can be -INF
 double f_cond_upper;            ///< conditional upper bound, can be +INF
 
 boost::multi_array< ColVariable , 2 > xi; ///< the observation variables
 std::vector< ColVariable > zeta;   ///< the satellite activation variables
 std::vector< ColVariable > activation; ///< the observation variables
 
 std::vector< FRowConstraint > orbitSelection; /// the satellite activation constraint

 boost::multi_array< FRowConstraint , 2 > activationSat_cnst; ///< the activation constraints 
 boost::multi_array< FRowConstraint , 2 > obs2_cnst; /// the (big-M linearized) observation constraints for CoverageSatLat 
 boost::multi_array< FRowConstraint , 2 > obs4_cnst; /// the (big-M linearized) observation constraints for CoverageSatLong

 FRealObjective c;               ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// register SatelliteBlock methods into the method factories
/** Although in general private methods should not be commented, this one is
 * because it does the registration of the following SatelliteBlock methods*/

 void guts_of_destructor( void );

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );


/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert SatelliteBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( SatelliteBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SatelliteBlockMod -----------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a SatelliteBlock
/** Derived class from Modification to describe modifications to a SatelliteBlock.
 *  This is actually "sort of abstract", since it does not say exactly what
 *  is changed, this being demanded to derived classes (which do this in
 *  different ways). Note that it is derived from Modification rather than,
 *  say, BlockMod (which has the same structure) because this is a class of
 *  "physical Modification". This means that a SatelliteBlockMod refers to changes
 *  in the "physical representation" of the SatelliteBlock; the corresponding
 *  changes in the "abstract representation" of the SatelliteBlock are dealt with
 *  by means of "abstract Modification", i.e., derived classes from
 *  AModification (as is BlockMod, which is why SatelliteBlockMod is not derived
 *  from BlockMod). */

class SatelliteBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of SatelliteBlockMod

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SatelliteBlock and the type

 SatelliteBlockMod( SatelliteBlock * fblock , int type )
  : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SatelliteBlockMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// returns the [DCR]Block to which the SatelliteBlockMod refers

 Block * get_Block( void ) const override  { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of modification

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SatelliteBlockMod

 void print( std::ostream &output ) const override {

  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 SatelliteBlock *f_Block;
               ///< pointer to the SatelliteBlock to which the SatelliteBlockMod refers

 int f_type;   ///< type of Modification

/*--------------------------------------------------------------------------*/

 };  // end( class( SatelliteBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SatelliteBlockRngdMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SatelliteBlockMod for "ranged" modifications
/** Derived class from SatelliteBlockMod to describe "ranged"
 * modifications to a SatelliteBlock.
 */

class SatelliteBlockRngdMod : public SatelliteBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SatelliteBlock, the type, and the range

 SatelliteBlockRngdMod( SatelliteBlock * fblock , int type , Block::Range rng )
  : SatelliteBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SatelliteBlockRngdMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return( f_rng ); }
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SatelliteBlockRngdMod

 void print( std::ostream &output ) const override {
  SatelliteBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;     ///< the range

/*--------------------------------------------------------------------------*/

 };  // end( class( SatelliteBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SatelliteBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SatelliteBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  SatelliteBlock. 
 */

class SatelliteBlockSbstMod : public SatelliteBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:


/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the SatelliteBlock, the type, and the subset
 /**< Constructor: takes the SatelliteBlock, the type, and the subset. As the the
  * && tells, nms is "consumed" by the constructor and its resources become
  * property of the SatelliteBlockSbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 SatelliteBlockSbstMod( SatelliteBlock * fblock , int type , Block::Subset && nms )
  : SatelliteBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SatelliteBlockSbstMod() = default;  ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SatelliteBlockSbstMod

 void print( std::ostream &output ) const override {
  SatelliteBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms;   ///< the subset

/*--------------------------------------------------------------------------*/

 };  // end( class( SatelliteBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SatelliteSolution -----------------------*/
/*--------------------------------------------------------------------------*/

class SatelliteSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*------------------------------- FRIENDS ----------------------------------*/

friend SatelliteBlock;  ///< make SatelliteBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING SatelliteSolution ----------*/

  explicit SatelliteSolution( void ) { }  /// constructor, it has nothing to do

  void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~SatelliteSolution() = default;  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A SatelliteSolution ------*/

  void read( const Block * block ) override final;

  void write( Block * block ) override final;

  void serialize( netCDF::NcGroup & group ) const override final;

  SatelliteSolution * scale( double factor ) const override final;

  void sum( const Solution * solution , double multiplier ) override final;

  SatelliteSolution * clone( bool empty = false ) const override final;
  
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

//protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override final {
   //output << "SatelliteSolution";
 }
  
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

//private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

SatelliteBlock::Vec_FNumber v_zeta;   ///< the arc flows
  
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( SatelliteSolution ) )

/** @} end( group( SatelliteBlock_CLASSES ) ) --------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SatelliteBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File SatelliteBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
