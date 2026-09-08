/*--------------------------------------------------------------------------*/
/*-------------------- File MultiTargetBlockv2.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiTargetBlockv2, which implements
 * the Block concept [see Block.h] for the solution of MultiTarget observavility problem.
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

#ifndef __MultiTargetBlockv2
#define __MultiTargetBlockv2 /* self-identification: #endif at the end of the file */

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
namespace SMSpp_di_unipi_it {
class MultiTargetBlockv2; // forward declaration of MultiTargetBlockv2

class MultiTargetSolution; // forward declaration of MultiTargetSolution

/*--------------------------------------------------------------------------*/
/*-------------------- MultiTargetBlockv2-RELATED TYPES --------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlockv2_TYPES MultiTargetBlockv2-related types
 *  @{ */

using p_MultiTargetBlockv2 =
 MultiTargetBlockv2 *; ///< a pointer to MultiTargetBlockv2

using Vec_MultiTargetBlockv2 = std::vector< p_MultiTargetBlockv2 >;
///< a vector of pointers to MultiTargetBlockv2

using Vec_MultiTargetBlockv2_it = Vec_MultiTargetBlockv2::iterator;
///< iterator for a Vec_MultiTargetBlockv2

using c_Vec_MultiTargetBlockv2 = const Vec_MultiTargetBlockv2;
///< a const vector of pointers to MultiTargetBlockv2

using c_Vec_MultiTargetBlockv2_it = c_Vec_MultiTargetBlockv2::iterator;
///< iterator for a c_Vec_MultiTargetBlockv2

/** @}  end( group( MultiTargetBlockv2_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlockv2_CLASSES Classes in MultiTargetBlockv2.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MultiTargetBlockv2 ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the "flat", non-decomposed
/// multi-target constellation-design problem
/** MultiTargetBlockv2 solves exactly the same problem as MultiTargetBlock
 * [see MultiTargetBlock.h], namely minimizing the average worst-case
 * revisit time over all targets subject to a common satellite orbit
 * configuration, but as a *single monolithic Block* rather than as one
 * MultiTargetBlock containing one SingleTargetBlock per target linked by
 * "duplicate" consistency constraints [see SingleTargetBlock.h for the
 * detailed mathematical model, which applies verbatim here]. Concretely:
 *
 * - the satellite orbit-selection Variable activation[ j ][ c ] (one per
 *   satellite j and candidate orbit c) is a *single* shared copy, since
 *   there is no decomposition by target requiring it to be duplicated and
 *   equated across per-target copies;
 *
 * - every other Variable and Constraint of SingleTargetBlock.h that used
 *   to live once per SingleTargetBlock (theta, Deltat, Deltat_k1/k2, zeta,
 *   b1/b2, d1/d2, h, xi, and all constraints built out of them) instead
 *   gains an explicit target dimension tgt, so that e.g. what was
 *   SingleTargetBlock::xi[ i ][ j ] (satellite i, time j, for one target)
 *   becomes MultiTargetBlockv2::xi[ i ][ j ][ m ] (satellite i, time j,
 *   target m);
 *
 * - consequently no "duplicate" constraints are needed at all: sharing is
 *   automatic since activation[][] is not replicated in the first place.
 *
 * Being avoided the overhead (and, for some Solver, the added complexity)
 * of nested sub-Block, MultiTargetBlockv2 is meant to be solved directly
 * by a general-purpose MILP Solver rather than through the decomposition
 * approach that MultiTargetBlock (and, similarly, ConstellationBlock with
 * its own SatelliteBlock) is designed for. Unlike its sibling classes,
 * MultiTargetBlockv2 also defines a full set of "physical Modification"
 * classes (MultiTargetBlockv2Mod and its Rngd/Sbst specializations, see
 * below), mirroring the ones of MCFBlock; however, as of now
 * guts_of_add_Modification() never actually issues one, so this machinery
 * is currently unused. */

class MultiTargetBlockv2 : public Block
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
 * MultiTargetBlockv2 defines three main public types:
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
 * However, while using a MultiTargetBlockv2 as a part of some larger problem, it may
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

 typedef double FNumber; ///< type of arc flow / deficit
 typedef const FNumber c_FNumber; ///< a read-only FNumber

 typedef std::vector< FNumber > Vec_FNumber; ///< a vector of FNumber
 typedef const Vec_FNumber c_Vec_FNumber; ///< a const vector of FNumber

 typedef Vec_FNumber::iterator Vec_FNumber_it; ///< iterator in Vec_FNumber
 typedef Vec_FNumber::const_iterator c_Vec_FNumber_it;
 ///< const iterator in Vec_FNumber

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double CNumber; ///< type of arc cost / potential
 typedef const CNumber c_CNumber; ///< a read-only CNumber

 typedef std::vector< CNumber > Vec_CNumber; ///< a vector of CNumber
 typedef const Vec_CNumber c_Vec_CNumber; ///< a const vector of CNumber

 typedef Vec_CNumber::iterator Vec_CNumber_it; ///< iterator in Vec_CNumber
 typedef Vec_CNumber::const_iterator c_Vec_CNumber_it;
 ///< const iterator in Vec_CNumber

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double FONumber;
 /**< type of the objective function: has to hold sums of products of
    FNumber(s) by CNumber(s) */

 typedef const FONumber c_FONumber; ///< a read-only FONumber

 typedef std::vector< FONumber > Vec_FONumber; ///< a vector of FONumber
 typedef const Vec_FONumber c_Vec_FONumber; ///< a const vector of FONumber

 /** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend MultiTargetSolution; ///< make MultiTargetSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructor and Destructor
 *  @{ */

 /// constructor of MultiTargetBlockv2, taking a pointer to the father (generic) Block
 /** Constructor of MultiTargetBlockv2. It accepts a pointer to the father Block, which
  * can be of any type, defaulting to nullptr so that this can also be used as
  * the void constructor. */

 explicit MultiTargetBlockv2( Block * father = nullptr )
     : Block( father ), AR1( 0 ), AR2( 0 ), AR3( 0 )
 {
 }


/*--------------------------------------------------------------------------*/
 /// destructor of MultiTargetBlockv2: deletes the abstract representation, if any

 virtual ~MultiTargetBlockv2() { guts_of_destructor(); }

 /** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Other initializations
 *  @{ */
 /*
  * Like load( std::istream & ), if there is any Solver attached to this
  * MultiTargetBlockv2 then a NBModification (the "nuclear option") is issued. */

 void load( const std::string & input , char frmt = 0 ) override;

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the DCR
 /** Method that generates the abstract Variable of the MultiTarget. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the MultiTarget
 /** Method that generates the abstract constraint of the MultiTarget. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the MultiTarget
 /** Method that generates the objective of the MultiTarget. */

 void generate_objective( Configuration * objc = nullptr ) override;

 //void generate_dynamic_constraints( Configuration *stcc = nullptr ) override;

 /** @} ---------------------------------------------------------------------*/
/*--------- Methods for reading the data of the MultiTargetBlockv2 ---------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for reading the data of the MultiTargetBlockv2
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override
 {
  return ( Objective::eMin );
 }

 /** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for checking the Block
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the MultiTargetBlockv2 is approximately feasible. This clearly
  * requires the Variable of the MultiTargetBlockv2 to have been defined, i.e., that
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

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) optimal
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the MultiTargetBlockv2 is approximately optimal, which means that
  * it is approximately feasible, that the dual solution encoded in the
  * current value of the dual multipliers of both the flow conservation and
  * bound constraints is approximately feasible, and that the two
  * approximately satisfies the Complementary Slackness Conditions. This
  * clearly requires that both the Variable and the Constraint of the
  * MultiTargetBlockv2 to have been defined, i.e., that generate_abstract_variables()
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

 bool is_optimal( bool useabstract = false ,
                  Configuration * optc = nullptr ) override;

 /** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for handling Solution
 *  @{ */

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// returns a MultiTargetSolution representing the current solution of this MultiTargetBlockv2

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

 /// returns the current value of Deltat[ k ], the revisit time of target k

 FNumber get_zeta( Index k ) const { return( Deltat[ k ].get_value() ); }

 /// sets the values of the Deltat ColVariable-s, one per target, in the
 /// range [ rng.first , rng.second ) from fstrt

 void set_zeta( c_Vec_FNumber_it fstrt ,
                Range rng = Range( 0 , Inf< Index >() ) );

 /** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for handling Modification
 *  @{ */

 /// returns true if there is any Solver "listening to this MultiTargetBlockv2"
 /** Returns true if there is any Solver "listening to this MultiTargetBlockv2", or if
  * the MultiTargetBlockv2 has to "listen" anyway because the "abstract" representation
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
 /// adding a new Modification to the MultiTargetBlockv2
 /** Method for handling Modification.
  *
  * The version of MultiTargetBlockv2 has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the MultiTargetBlockv2, and "translate"
  * them into both changes of the actual data structures and corresponding
  * "physical Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the MultiTargetBlockv2, and this must not be
  * repeated. In particular, this would be an issue if the Modification would
  * be [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * MultiTargetBlockv2), which would mean that the "physical Modification" would be
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
  * Any other Modification reaching the MultiTargetBlockv2 will lead to exception
  * being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

 /** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE MultiTargetBlockv2 ----------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for printing & saving the MultiTargetBlockv2
 *  @{ */

 /// print the MultiTargetBlockv2 on an ostream with the given verbosity
 /** Protected method to print information about the MultiTargetBlockv2; with the
  * "complete" level ('C') it outputs the MultiTargetBlockv2 in DIMACS format. */

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

 unsigned char AR1, AR2, AR3; ///< bit-wise coded: what abstract is there
 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variable have been constructed
 static constexpr unsigned char HasObj = 2;
 ///< first bit of AR == 1 if the Objective have been constructed
 static constexpr unsigned char HasCnst = 3;
 ///< first bit of AR == 1 if the Constraint have been constructed

 Index n; ///< the number of satellites
 Index t; ///< the total number of time step
 FNumber OrbitSet; ///< the number of configurations

 Index tgt; ///< the total number of time step
 FNumber targets;

 FNumber indexOrbit;
 FNumber dt; ///< the time discretization step
 FNumber T; ///< the simulation horizon
 FNumber alphaHalf; ///< the alpha half MultiTarget parameter

 FNumber altitudeVal; ///< vector of MultiTarget altitude values
 FNumber thetaValFinal; ///< vector of MultiTarget theta values

 FNumber satellites; ///< the number of satellites in the MultiTarget
 FNumber time_step; ///< the time discretization step
 FNumber horizon; ///< the simulation horizon [sec]

 boost::multi_array< double, 3 >
  CoverageSatLat; ///< the matrix of pre-computed CoverageSatLat
 boost::multi_array< double, 3 >
  CoverageSatLong; ///< the matrix of pre-computed CoverageSatLong

 double f_cond_lower; ///< conditional lower bound, can be -INF
 double f_cond_upper; ///< conditional upper bound, can be +INF

 std::vector< ColVariable > theta; ///< the satellite theta

 std::vector< ColVariable > Deltat; ///< the MultiTarget revisit time
 boost::multi_array< ColVariable, 2 >
  Deltat_k1; ///< the MultiTarget revisit time (1)
 boost::multi_array< ColVariable, 2 >
  Deltat_k2; ///< the MultiTarget revisit time (2)
 boost::multi_array< ColVariable, 2 >
  zeta; ///< the MultiTarget activation variables
 boost::multi_array< ColVariable, 2 >
  b1; ///< the MultiTarget activation variables
 boost::multi_array< ColVariable, 2 >
  b2; ///< the MultiTarget activation variables

 boost::multi_array< ColVariable, 2 > Deltat_k1A;
 boost::multi_array< ColVariable, 2 > Deltat_k2A;
 boost::multi_array< ColVariable, 2 > Deltat_k2AUX;

 boost::multi_array< ColVariable, 3 > xi; ///< the observation variables
 boost::multi_array< ColVariable, 2 > h; ///< the two-observation variables
 boost::multi_array< ColVariable, 2 >
  activation; ///< the observation variables

 boost::multi_array< ColVariable, 2 > d1; ///< the d1 variables
 boost::multi_array< ColVariable, 2 > d2; ///< the d2 variables

 std::vector< FRowConstraint >
  orbitSelection; /// the MultiTarget activation constraint
 std::vector< FRowConstraint >
  theta_UB; /// the MultiTarget activation constraint

 //boost::multi_array< FRowConstraint , 3 >  duplicate_pi; /// the duplicate_pi constraints

 boost::multi_array< FRowConstraint, 2 >
  observation1; /// the observation1 constraints

 boost::multi_array< FRowConstraint, 2 >
  Deltat_max1; /// the Deltat_max1 constraints
 boost::multi_array< FRowConstraint, 2 >
  Deltat_max11; /// the Deltat_max1 constraints
 boost::multi_array< FRowConstraint, 2 >
  Deltat_max2; /// the Deltat_max2 constraints
 boost::multi_array< FRowConstraint, 2 >
  Deltat_max22; /// the Deltat_max1 constraints
 boost::multi_array< FRowConstraint, 2 >
  Deltat_max3; /// the Deltat_max3 constraints

 boost::multi_array< FRowConstraint, 2 >
  Deltat_min_k1_1; ///< the Deltat_min_k1 constraints (1)
 boost::multi_array< FRowConstraint, 2 >
  Deltat_min_k1_2; ///< the Deltat_min_k1 constraints (2)
 boost::multi_array< FRowConstraint, 2 >
  d1_cnst; /// the d1 activation constraint

 boost::multi_array< FRowConstraint, 2 >
  Deltat_min_k1_11; ///< the Deltat_min_k1 constraints (1)
 boost::multi_array< FRowConstraint, 2 > Deltat_min_k2_11;

 boost::multi_array< FRowConstraint, 2 >
  Deltat_min_k2_1; ///< the Deltat_min_k2 constraints (1)
 boost::multi_array< FRowConstraint, 2 >
  Deltat_min_k2_2; ///< the Deltat_min_k2 constraints (2)
 boost::multi_array< FRowConstraint, 2 >
  d2_cnst; /// the d2 activation constraint
 boost::multi_array< FRowConstraint, 2 >
  d2a_cnst; /// the d2 activation constraint

 boost::multi_array< FRowConstraint, 2 > d1A_cnst;
 boost::multi_array< FRowConstraint, 2 > d2A_cnst;

 boost::multi_array< FRowConstraint, 3 >
  activationSat_cnst; ///< the observation constraints
 boost::multi_array< FRowConstraint, 2 >
  activationSat1_cnst; ///< the observation constraints

 boost::multi_array< FRowConstraint, 2 >
  h_cnst_1; ///< the linearization of the product of two zeta's (1)
 boost::multi_array< FRowConstraint, 2 >
  h_cnst_2; ///< the linearization of the product of two zeta's (2)
 boost::multi_array< FRowConstraint, 2 >
  h_cnst_3; ///< the linearization of the product of two zeta's (3)

 boost::multi_array< FRowConstraint, 3 >
  obs1_cnst; /// the linearized observation constraints via big-M for CoverageSatLat (1)
 boost::multi_array< FRowConstraint, 3 >
  obs2_cnst; /// the linearized observation constraints via big-M for CoverageSatLat (2)
 boost::multi_array< FRowConstraint, 3 >
  obs3_cnst; /// the linearized observation constraints via big-M for CoverageSatLong (1)
 boost::multi_array< FRowConstraint, 3 >
  obs4_cnst; /// the linearized observation constraints via big-M for CoverageSatLong (2)

 std::vector< FRowConstraint > obs_cnst;
 std::vector< FRowConstraint > obs_cnst_h;
 boost::multi_array< FRowConstraint, 2 > obs_cnst_xi;

 std::vector< FRowConstraint > Deltat_max_dt;
 std::vector< FRowConstraint > Deltat_max_dt1;

 FRealObjective c; ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:
/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
 /// register MultiTargetBlockv2 methods into the method factories
 /** Although in general private methods should not be commented, this one is
 * because it does the registration of the following MultiTargetBlockv2 methods*/

 void guts_of_destructor( void );

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );


/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h; // insert MultiTargetBlockv2 in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}; // end( class( MultiTargetBlockv2 ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS MultiTargetBlockv2Mod -----------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a MultiTargetBlockv2
/** Derived class from Modification to describe modifications to a MultiTargetBlockv2.
 *  This is actually "sort of abstract", since it does not say exactly what
 *  is changed, this being demanded to derived classes (which do this in
 *  different ways). Note that it is derived from Modification rather than,
 *  say, BlockMod (which has the same structure) because this is a class of
 *  "physical Modification". This means that a MultiTargetBlockv2Mod refers to changes
 *  in the "physical representation" of the MultiTargetBlockv2; the corresponding
 *  changes in the "abstract representation" of the MultiTargetBlockv2 are dealt with
 *  by means of "abstract Modification", i.e., derived classes from
 *  AModification (as is BlockMod, which is why MultiTargetBlockv2Mod is not derived
 *  from BlockMod).
 *
 * NOTE: as of now, MultiTargetBlockv2::guts_of_add_Modification() never
 * actually constructs any MultiTargetBlockv2Mod (or its Rngd/Sbst
 * specializations below), so this class is currently unused. */

class MultiTargetBlockv2Mod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

public:
/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of MultiTargetBlockv2Mod

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the MultiTargetBlockv2 and the type

 MultiTargetBlockv2Mod( MultiTargetBlockv2 * fblock , int type )
     : f_Block( fblock ), f_type( type )
 {
 }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~MultiTargetBlockv2Mod() = default; ///< destructor, does nothing

/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/

 /// returns the [DCR]Block to which the MultiTargetBlockv2Mod refers

 Block * get_Block( void ) const override { return ( f_Block ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of modification

 int type( void ) const { return ( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

protected:
/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the MultiTargetBlockv2Mod

 void print( std::ostream & output ) const override {}

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 MultiTargetBlockv2 * f_Block;
 ///< pointer to the MultiTargetBlockv2 to which the MultiTargetBlockv2Mod refers

 int f_type; ///< type of Modification

/*--------------------------------------------------------------------------*/

}; // end( class( MultiTargetBlockv2Mod ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS MultiTargetBlockv2RngdMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from MultiTargetBlockv2Mod for "ranged" modifications
/** Derived class from MultiTargetBlockv2Mod to describe "ranged"
 * modifications to a MultiTargetBlockv2, i.e., modifications that apply to an interval
 * of either arcs or nodes. */

class MultiTargetBlockv2RngdMod : public MultiTargetBlockv2Mod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

public:
/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the MultiTargetBlockv2, the type, and the range

 MultiTargetBlockv2RngdMod( MultiTargetBlockv2 * fblock , int type ,
                            Block::Range rng )
     : MultiTargetBlockv2Mod( fblock , type ), f_rng( rng )
 {
 }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~MultiTargetBlockv2RngdMod() = default; ///< destructor, does nothing

/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return ( f_rng ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

protected:
/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the MultiTargetBlockv2RngdMod

 void print( std::ostream & output ) const override
 {
  MultiTargetBlockv2Mod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng; ///< the range

/*--------------------------------------------------------------------------*/

}; // end( class( MultiTargetBlockv2RngdMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS MultiTargetBlockv2SbstMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from MultiTargetBlockv2Mod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  MultiTargetBlockv2, i.e., modifications that apply to an arbitrary subset of either
 * the arcs or the nodes. */

class MultiTargetBlockv2SbstMod : public MultiTargetBlockv2Mod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

public:
/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the MultiTargetBlockv2, the type, and the subset
 /**< Constructor: takes the MultiTargetBlockv2, the type, and the subset. As the the
  * && tells, nms is "consumed" by the constructor and its resources become
  * property of the MultiTargetBlockv2SbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 MultiTargetBlockv2SbstMod( MultiTargetBlockv2 * fblock , int type ,
                            Block::Subset && nms )
     : MultiTargetBlockv2Mod( fblock , type ), f_nms( std::move( nms ) )
 {
 }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~MultiTargetBlockv2SbstMod() = default; ///< destructor, does nothing

/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return ( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

protected:
/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the MultiTargetBlockv2SbstMod

 void print( std::ostream & output ) const override
 {
  MultiTargetBlockv2Mod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms; ///< the subset

/*--------------------------------------------------------------------------*/

}; // end( class( MultiTargetBlockv2SbstMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS MultiTargetSolution ------------------------*/
/*--------------------------------------------------------------------------*/

class MultiTargetSolution : public Solution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*------------------------------- FRIENDS ----------------------------------*/

 friend MultiTargetBlockv2; ///< make MultiTargetBlockv2 friend

/*------------ CONSTRUCTING AND DESTRUCTING MultiTargetSolution ------------*/

 explicit MultiTargetSolution( void ) {} /// constructor, it has nothing to do

 void deserialize( const netCDF::NcGroup & group ) override final;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~MultiTargetSolution() = default; ///< destructor: it is virtual, and empty

/*-------- METHODS DESCRIBING THE BEHAVIOR OF A MultiTargetSolution --------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

 void serialize( netCDF::NcGroup & group ) const override final;

 MultiTargetSolution * scale( double factor ) const override final;

 void sum( const Solution * solution , double multiplier ) override final;

 MultiTargetSolution * clone( bool empty = false ) const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 //protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override final
 {
  //output << "MultiTargetSolution";
 }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 //private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 MultiTargetBlockv2::Vec_FNumber v_zeta; ///< the arc flows

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

}; // end( class( MultiTargetSolution ) )

/** @} end( group( MultiTargetBlockv2_CLASSES ) ) --------------------------*/
/*--------------------------------------------------------------------------*/

}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* MultiTargetBlockv2.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File MultiTargetBlockv2.h ------------------------*/
/*--------------------------------------------------------------------------*/
