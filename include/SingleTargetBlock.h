/*--------------------------------------------------------------------------*/
/*-------------------- File SingleTargetBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SingleTargetBlock, which implements
 * the Block concept [see Block.h] for the solution of SingleTarget observavility problem.
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

#ifndef __SingleTargetBlock
 #define __SingleTargetBlock  /* self-identification: #endif at the end of the file */

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
 class SingleTargetBlock;     // forward declaration of SingleTargetBlock

 class SingleTargetSolution;  // forward declaration of SingleTargetSolution

/*--------------------------------------------------------------------------*/
/*----------------------- SingleTargetBlock-RELATED TYPES ------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleTargetBlock_TYPES SingleTargetBlock-related types
 *  @{ */

 using p_SingleTargetBlock = SingleTargetBlock *;  ///< a pointer to SingleTargetBlock

 using Vec_SingleTargetBlock = std::vector< p_SingleTargetBlock>;
 ///< a vector of pointers to SingleTargetBlock

 using Vec_SingleTargetBlock_it = Vec_SingleTargetBlock::iterator;
 ///< iterator for a Vec_SingleTargetBlock

 using c_Vec_SingleTargetBlock = const Vec_SingleTargetBlock;
 ///< a const vector of pointers to SingleTargetBlock

 using c_Vec_SingleTargetBlock_it = c_Vec_SingleTargetBlock::iterator;
 ///< iterator for a c_Vec_SingleTargetBlock

/** @}  end( group( SingleTargetBlock_TYPES ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SingleTargetBlock_CLASSES Classes in SingleTargetBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleTargetBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple SingleTargetBlock concept.
/*
* Let T and be the time horizon, we indicate with [T] := \{ 1,2,...,T \} 
* the set of time stamps. dt is the time-step for time discretization. 
* [C] is the set of possible orbital configurations.
*
* We consider two sets of non-negative continuous decision variables:
*
* (1) \Delta_t is the maximum revisit time for the current target
* (2) \Delta_t^1[ k ] and \Delta_t^2[ k ] are are auxiliary variables to 
* define the revisit times before time stamp k \in [T] for the current target
*
* We consider four sets of binary decision variables:
*
* (1) \pi[ c ] \in \{ 0 , 1 \}, c \in [C], indicating which orbit configuration 
* c \in [ C ] is selected for the satellite, 
* (2) \xi[ t ][ m ] \in \{ 0 , 1 \} indicating whether the satellite i observes 
* the curent target at time t, and
* (3) zeta[ t ] \in \{ 0 , 1 \} indicating whether the constellation observes  
* the current target
* (4) h[ k ][ j ] \in \{ 0 , 1 \} indicating whether the constellation observes  
* the current target in two time stamps k and j
* 
* The model we are going to solve reads as follows.
*
* \f[
*   min \Delta_t                                                    (1)
* \f]
* \f[
*   min \Delta_t^1[ k ] = \min_{j \in [T] : j \leq k} 
*   \{ dt * (k-j) * h[ k ][ j ], + T/2 * ( 1 - h[ k ][ j ] )\}, 
*   \forall k \in [T]                                               (2)
* \f]
*   min \Delta_t^2[ k ] = \min_{j \in [T] : j \leq k} 
*   \{ T - dt * (k-j) * h[ k ][ j ], + T/2 * ( 1 - h[ k ][ j ] )\}, 
*   \forall k \in [T]                                               (3)
* \f]
* \f[
*   min \Delta_t >= \max \{ \Delta_t^1[ k ], \Delta_t^2[ k ] \}     (4)
* \f]
* \f[
* \xi[ t ][ m ]=0 \Longrightarrow \max \{ \sum_{c\in [C]} \pi_{c}\, 
\Delta lat[ c ][ t ][ m ]}, \sum_{c\in [C]} \pi_{cc} \,\Delta long[ c ][ t ][ m ] \} 
\ge \theta^{\max},       \forall t\in T(dt), \forall m \in \mathcal{X}     (1)
* \f]
*
* The objective function (1) minimizes the maximum revisit time for the current 
* target Constraints (2) and (3) define \Delta_t^1[ k ] and \Delta_t^2[ k ],
* respectively, as the revisit times before time stamp k \in [T] for the current 
* target (these constraints are opportunely linearized when defining the constraint
* of the Block). Then, the maximum revisit time is the maximum between \Delta_t^1[ k ]
* and \Delta_t^2[ k ], see constraint (4).
*/

class SingleTargetBlock : public Block
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
 * SingleTargetBlock defines three main public types:
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
 * However, while using a SingleTargetBlock as a part of some larger problem, it may
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

 friend SingleTargetSolution;  ///< make SingleTargetSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of SingleTargetBlock, taking a pointer to the father (generic) Block
 /** Constructor of SingleTargetBlock. It accepts a pointer to the father Block, which
  * can be of any type, defaulting to nullptr so that this can also be used as
  * the void constructor. */

 explicit SingleTargetBlock( Block *father = nullptr )
  : Block( father ) , AR1(0), AR2(0), AR3(0) { }
               

/*--------------------------------------------------------------------------*/
 /// destructor of SingleTargetBlock: deletes the abstract representation, if any

 virtual ~SingleTargetBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
  /*
  * Like load( std::istream & ), if there is any Solver attached to this
  * SingleTargetBlock then a NBModification (the "nuclear option") is issued. */

 void load( FNumber n , FNumber dt , FNumber T , FNumber altValues , FNumber thetaValues ,
              FNumber indexOrbit, FNumber aHalf ,
              boost::multi_array< double , 2 > CoverageLat , boost::multi_array< double , 2 > CoverageLong );

 void load( std::istream &input , char frmt = 0 ) override;

 [[nodiscard]] ColVariable * i2p_pi( Index iii , Index jjj ) const {
  return( const_cast< ColVariable * >( &activation[ iii ][ jjj ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the SingleTargetBlock
 /** Method that generates the abstract Variable of the SingleTarget. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the SingleTarget
 /** Method that generates the abstract constraint of the SingleTarget. */
 
 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the SingleTarget
 /** Method that generates the objective of the SingleTarget. */

 void generate_objective( Configuration *objc = nullptr ) override;

 //void generate_dynamic_constraints( Configuration *stcc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SingleTargetBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the SingleTargetBlock
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
  * (x) Variable of the SingleTargetBlock is approximately feasible. This clearly
  * requires the Variable of the SingleTargetBlock to have been defined, i.e., that
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
  * (x) Variable of the SingleTargetBlock is approximately optimal, which means that
  * it is approximately feasible, that the dual solution encoded in the
  * current value of the dual multipliers of both the flow conservation and
  * bound constraints is approximately feasible, and that the two
  * approximately satisfies the Complementary Slackness Conditions. This
  * clearly requires that both the Variable and the Constraint of the
  * SingleTargetBlock to have been defined, i.e., that generate_abstract_variables()
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
  
 /// returns a SingleTargetSolution representing the current solution of this SingleTargetBlock

 Solution * get_Solution( Configuration *solc = nullptr ,
 			  bool emptys = true ) override;


/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// returns true if there is any Solver "listening to this SingleTargetBlock"
 /** Returns true if there is any Solver "listening to this SingleTargetBlock", or if
  * the SingleTargetBlock has to "listen" anyway because the "abstract" representation
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
 /// adding a new Modification to the SingleTargetBlock
 /** Method for handling Modification.
  *
  * The version of SingleTargetBlock has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the SingleTargetBlock, and "translate"
  * them into both changes of the actual data structures and corresponding
  * "physical Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the SingleTargetBlock, and this must not be
  * repeated. In particular, this would be an issue if the Modification would
  * be [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * SingleTargetBlock), which would mean that the "physical Modification" would be
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
  * Any other Modification reaching the SingleTargetBlock will lead to exception
  * being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE SingleTargetBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SingleTargetBlock
 *  @{ */

 /// print the SingleTargetBlock on an ostream with the given verbosity
 /** Protected method to print information about the SingleTargetBlock; with the
  * "complete" level ('C') it outputs the SingleTargetBlock in DIMACS format. */

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

 Index n;                        ///< the number of satellites
 Index t;                        ///< the total number of time step
 FNumber OrbitSet;               ///< the number of configurations

 FNumber indexOrbit;
 FNumber dt;                     ///< the time discretization step
 FNumber T;                      ///< the simulation horizon
 FNumber alphaHalf;              ///< the alpha half SingleTarget parameter

 FNumber altitudeVal;        ///< vector of satellite altitude values
 FNumber thetaVal;           ///< vector of satellite theta values

 boost::multi_array< double , 2 > CoverageSatLat;  ///< the matrix of pre-computed CoverageSatLat
 boost::multi_array< double , 2 > CoverageSatLong; ///< the matrix of pre-computed CoverageSatLong

 double f_cond_lower;            ///< conditional lower bound, can be -INF
 double f_cond_upper;            ///< conditional upper bound, can be +INF
 
 std::vector< ColVariable > Deltat;   ///< the SingleTarget revisit time
 std::vector< ColVariable > Deltat_k1;   ///< the SingleTarget revisit time (1)
 std::vector< ColVariable > Deltat_k2;   ///< the SingleTarget revisit time (2)
 std::vector< ColVariable > zeta;   ///< the SingleTarget activation variables
 std::vector< ColVariable > b1;   ///< the SingleTarget activation variables
 std::vector< ColVariable > b2;   ///< the SingleTarget activation variables

 std::vector< ColVariable > Deltat_k1A;
 std::vector< ColVariable > Deltat_k2A;
 std::vector< ColVariable > Deltat_k2AUX;

 boost::multi_array< ColVariable , 2 > xi; ///< the observation variables
 std::vector< ColVariable > h; ///< the two-observation variables
 boost::multi_array< ColVariable , 2 > activation; ///< the observation variables

 std::vector< ColVariable > d1; ///< the d1 variables
 std::vector< ColVariable > d2; ///< the d2 variables
  
 std::vector< FRowConstraint > orbitSelection; /// the SingleTarget activation constraint

 std::vector< FRowConstraint > Deltat_max1; /// the Deltat_max1 constraints
 std::vector< FRowConstraint > Deltat_max11; /// the Deltat_max1 constraints
 std::vector< FRowConstraint > Deltat_max2; /// the Deltat_max2 constraints
 std::vector< FRowConstraint > Deltat_max22; /// the Deltat_max1 constraints
 std::vector< FRowConstraint > Deltat_max3; /// the Deltat_max3 constraints

 std::vector< FRowConstraint > Deltat_min_k1_1; ///< the Deltat_min_k1 constraints (1)
 std::vector< FRowConstraint > Deltat_min_k1_2; ///< the Deltat_min_k1 constraints (2)
 std::vector< FRowConstraint > d1_cnst; /// the d1 activation constraint

 std::vector< FRowConstraint > Deltat_min_k1_11; ///< the Deltat_min_k1 constraints (1)
 std::vector< FRowConstraint > Deltat_min_k2_11;

 std::vector< FRowConstraint > Deltat_min_k2_1; ///< the Deltat_min_k2 constraints (1)
 std::vector< FRowConstraint > Deltat_min_k2_2; ///< the Deltat_min_k2 constraints (2)
 std::vector< FRowConstraint > d2_cnst; /// the d2 activation constraint
 std::vector< FRowConstraint > d2a_cnst; /// the d2 activation constraint
  
 std::vector< FRowConstraint > d1A_cnst;
 std::vector< FRowConstraint > d2A_cnst;
  
 boost::multi_array< FRowConstraint , 2 > activationSat_cnst; ///< the observation constraints
 std::vector< FRowConstraint > activationSat1_cnst; ///< the observation constraints

 std::vector< FRowConstraint > h_cnst_1; ///< the linearization of the product of two zeta's (1)
 std::vector< FRowConstraint > h_cnst_2; ///< the linearization of the product of two zeta's (2)
 std::vector< FRowConstraint > h_cnst_3; ///< the linearization of the product of two zeta's (3)

 boost::multi_array< FRowConstraint , 2 > obs1_cnst; /// the linearized observation constraints via big-M for CoverageSatLat (1)
 boost::multi_array< FRowConstraint , 2 > obs2_cnst; /// the linearized observation constraints via big-M for CoverageSatLat (2)
 boost::multi_array< FRowConstraint , 2 > obs3_cnst; /// the linearized observation constraints via big-M for CoverageSatLong (1)
 boost::multi_array< FRowConstraint , 2 > obs4_cnst; /// the linearized observation constraints via big-M for CoverageSatLong (2)

 std::vector< FRowConstraint > obs_cnst;  
 std::vector< FRowConstraint > obs_cnst_h;  
 std::vector< FRowConstraint > obs_cnst_xi;  

 std::vector< FRowConstraint > Deltat_max_dt;
 std::vector< FRowConstraint > Deltat_max_dt1;
  
 FRealObjective c;               ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// register SingleTargetBlock methods into the method factories
/** Although in general private methods should not be commented, this one is
 * because it does the registration of the following SingleTargetBlock methods*/

 void guts_of_destructor( void );

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );


/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert SingleTargetBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( SingleTargetBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleTargetBlockMod -----------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a SingleTargetBlock
/** Derived class from Modification to describe modifications to a SingleTargetBlock.
 *  This is actually "sort of abstract", since it does not say exactly what
 *  is changed, this being demanded to derived classes (which do this in
 *  different ways). Note that it is derived from Modification rather than,
 *  say, BlockMod (which has the same structure) because this is a class of
 *  "physical Modification". This means that a SingleTargetBlockMod refers to changes
 *  in the "physical representation" of the SingleTargetBlock; the corresponding
 *  changes in the "abstract representation" of the SingleTargetBlock are dealt with
 *  by means of "abstract Modification", i.e., derived classes from
 *  AModification (as is BlockMod, which is why SingleTargetBlockMod is not derived
 *  from BlockMod). */

class SingleTargetBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of SingleTargetBlockMod

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SingleTargetBlock and the type

 SingleTargetBlockMod( SingleTargetBlock * fblock , int type )
  : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleTargetBlockMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// returns the [DCR]Block to which the SingleTargetBlockMod refers

 Block * get_Block( void ) const override  { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of modification

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SingleTargetBlockMod

 void print( std::ostream &output ) const override {

  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 SingleTargetBlock *f_Block;
               ///< pointer to the SingleTargetBlock to which the SingleTargetBlockMod refers

 int f_type;   ///< type of Modification

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleTargetBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SingleTargetBlockRngdMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SingleTargetBlockMod for "ranged" modifications
/** Derived class from SingleTargetBlockMod to describe "ranged"
 * modifications to a SingleTargetBlock. 
 */

class SingleTargetBlockRngdMod : public SingleTargetBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the SingleTargetBlock, the type, and the range

 SingleTargetBlockRngdMod( SingleTargetBlock * fblock , int type , Block::Range rng )
  : SingleTargetBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleTargetBlockRngdMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return( f_rng ); }
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SingleTargetBlockRngdMod

 void print( std::ostream &output ) const override {
  SingleTargetBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;     ///< the range

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleTargetBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SingleTargetBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from SingleTargetBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  SingleTargetBlock. 
 */

class SingleTargetBlockSbstMod : public SingleTargetBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the SingleTargetBlock, the type, and the subset
 /**< Constructor: takes the SingleTargetBlock, the type, and the subset. As the the
  * && tells, nms is "consumed" by the constructor and its resources become
  * property of the SingleTargetBlockSbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 SingleTargetBlockSbstMod( SingleTargetBlock * fblock , int type , Block::Subset && nms )
  : SingleTargetBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~SingleTargetBlockSbstMod() = default;  ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the SingleTargetBlockSbstMod

 void print( std::ostream &output ) const override {
  SingleTargetBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms;   ///< the subset

/*--------------------------------------------------------------------------*/

 };  // end( class( SingleTargetBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SingleTargetSolution --------------------*/
/*--------------------------------------------------------------------------*/

class SingleTargetSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*------------------------------- FRIENDS ----------------------------------*/

friend SingleTargetBlock;  ///< make SingleTargetBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING SingleTargetSolution ----------------*/

  explicit SingleTargetSolution( void ) { }  /// constructor, it has nothing to do

  void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~SingleTargetSolution() = default;  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A SingleTargetSolution -----------*/

  void read( const Block * block ) override final;

  void write( Block * block ) override final;

  void serialize( netCDF::NcGroup & group ) const override final;

  SingleTargetSolution * scale( double factor ) const override final;

  void sum( const Solution * solution , double multiplier ) override final;

  SingleTargetSolution * clone( bool empty = false ) const override final;
  
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

//protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override final {
   //output << "SingleTargetSolution";
 }
  
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

//private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

SingleTargetBlock::Vec_FNumber v_zeta;   ///< the arc flows
  
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( SingleTargetSolution ) )

/** @} end( group( SingleTargetBlock_CLASSES ) ) --------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SingleTargetBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File SingleTargetBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
