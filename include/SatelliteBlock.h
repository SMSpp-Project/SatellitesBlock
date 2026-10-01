/*--------------------------------------------------------------------------*/
/*------------------------- File SatelliteBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SatelliteBlock, which implements the
 * Block concept [see Block.h] for a single satellite of the Satellite
 * Constellation Design Problem, and for the class SatelliteSolution
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

#ifndef __SatelliteBlock
 #define __SatelliteBlock
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
class SatelliteBlock; // forward declaration of SatelliteBlock

class SatelliteSolution; // forward declaration of SatelliteSolution

/*--------------------------------------------------------------------------*/
/*---------------------- SatelliteBlock-RELATED TYPES ----------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteBlock_TYPES SatelliteBlock-related types
 *  @{ */

using p_SatelliteBlock = SatelliteBlock *; ///< a pointer to SatelliteBlock

using Vec_SatelliteBlock = std::vector< p_SatelliteBlock >;
///< a vector of pointers to SatelliteBlock

using Vec_SatelliteBlock_it = Vec_SatelliteBlock::iterator;
///< iterator for a Vec_SatelliteBlock

using c_Vec_SatelliteBlock = const Vec_SatelliteBlock;
///< a const vector of pointers to SatelliteBlock

using c_Vec_SatelliteBlock_it = c_Vec_SatelliteBlock::iterator;
///< iterator for a c_Vec_SatelliteBlock

/** @} end( group( SatelliteBlock_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteBlock_CLASSES Classes in SatelliteBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SatelliteBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a single satellite of a ConstellationBlock
/** SatelliteBlock represents a single satellite of the Satellite
 * Constellation Design Problem [see ConstellationBlock.h]. Let T be the
 * time horizon, dt the time step, \f$ \mathcal{X} \f$ the set of the n
 * targets, and [C] the set of the candidate orbits of the satellite, for
 * each of which the latitude and longitude distances
 * \f$ \Delta lat[ c ][ m ][ t ] \f$ and \f$ \Delta long[ c ][ m ][ t ] \f$
 * between the ground track and every target m at every time stamp t are
 * known (they are computed by ConstellationBlock::load()). The Variable are
 *
 * - \f$ \pi[ c ] \in \{ 0 , 1 \} \f$, c \in [C], which is 1 if the
 *   satellite is placed on orbit c (activation);
 *
 * - \f$ \xi[ m ][ t ] \in \{ 0 , 1 \} \f$, which is 1 if the satellite
 *   observes target m at time stamp t (xi);
 *
 * - \f$ \zeta \in \{ 0 , 1 \} \f$, which is 1 if the satellite is active
 *   (zeta);
 *
 * - \f$ \theta \geq 0 \f$, the observability threshold of the satellite
 *   (thetaVar).
 *
 * The satellite observes target m at time t only if both distances, for
 * the chosen orbit, are within the threshold:
 *
 * \f[
 *  \xi[ m ][ t ] = 1 \Longrightarrow \max \Bigl\{ \sum_{ c \in [C] }
 *  \pi[ c ] \Delta lat[ c ][ m ][ t ] \, , \, \sum_{ c \in [C] } \pi[ c ]
 *  \Delta long[ c ][ m ][ t ] \Bigr\} \leq \theta                 \qquad (1)
 * \f]
 *
 * which is linearized with big-M constraints (obs2_cnst for the latitude,
 * obs4_cnst for the longitude), with M the largest distance over [C]. An
 * active satellite has exactly one orbit, an inactive one none, and it is
 * active iff it observes some target:
 *
 * \f[
 *  \sum_{ c \in [C] } \pi[ c ] = \zeta                            \qquad (2)
 * \f]
 * \f[
 *  \xi[ m ][ t ] \leq \zeta \;\; \forall m , t \; , \quad
 *  \zeta \leq \sum_{ m , t } \xi[ m ][ t ]                        \qquad (3)
 * \f]
 *
 * Finally, the threshold is \f$ \theta^{\max} \zeta / 20 \leq \theta \leq
 * \theta^{\max} \zeta \f$, where \f$ \theta^{\max} \f$ is the maximum
 * threshold given to load(). The objective is \f$ \zeta \f$, so that the
 * sum of the objectives of the SatelliteBlock of a ConstellationBlock is
 * the number of active satellites. */

class SatelliteBlock : public Block {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * SatelliteBlock defines three main public types:
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

 friend SatelliteSolution; ///< make SatelliteSolution friend

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of SatelliteBlock, taking a pointer to the father Block
 /** Constructor of SatelliteBlock. It accepts a pointer to the father
  * Block, which can be of any type, defaulting to nullptr so that this can
  * also be used as the void constructor. */

 explicit SatelliteBlock( Block * father = nullptr )
  : Block( father ) , AR1( 0 ) , AR2( 0 ) , AR3( 0 ) , n( 0 ) , t( 0 ) ,
    OrbitSet( 0 ) , dt( 0 ) , T( 0 ) , alphaHalf( 0 ) , altitudeVal( 0 ) ,
    thetaVal( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of SatelliteBlock: deletes the abstract representation

 virtual ~SatelliteBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance out of memory
 /** Loads the instance out of memory; the data is computed by
  * ConstellationBlock::load(), which reads them out of a file. The
  * parameters are:
  *
  * - num_targets, the number n of targets;
  *
  * - time_step, the time step dt [s];
  *
  * - horizon, the time horizon T [s], so that there are T / dt time
  *   stamps;
  *
  * - altValues, the altitude of the orbits of the satellite [m];
  *
  * - thetaValues, the maximum threshold \f$ \theta^{\max} \f$;
  *
  * - indOrbit, the number of candidate orbits [C];
  *
  * - aHalf, the half-aperture of the sensor of the satellite;
  *
  * - CoverageLat and CoverageLong, with CoverageLat[ m ][ t ][ c ] the
  *   latitude distance between the ground track of orbit c and target m at
  *   time stamp t, and the same for the longitude; only the first indOrbit
  *   entries of the last dimension are used;
  *
  * - periods_tgt, the number of revisit periods of each target.
  *
  * Any previous abstract representation is deleted. If there is any Solver
  * attached to this SatelliteBlock then a NBModification (the "nuclear
  * option") is issued. */

 void load( Index num_targets , FNumber time_step , FNumber horizon ,
            FNumber altValues , FNumber thetaValues , Index indOrbit ,
            FNumber aHalf ,
            const boost::multi_array< double , 3 > & CoverageLat ,
            const boost::multi_array< double , 3 > & CoverageLong ,
            const std::vector< double > & periods_tgt );

/*--------------------------------------------------------------------------*/
 /// loading a SatelliteBlock out of an istream is not supported
 /** A SatelliteBlock is only loaded by its ConstellationBlock, out of
  * memory: this always throws std::logic_error. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// loads the instance from the file with the given name
 /** The version of Block, which opens the file and calls
  * load( std::istream & ). */

 using Block::load;

/*--------------------------------------------------------------------------*/
 /// generate the abstract Variable of the SatelliteBlock
 /** Generates the Variable thetaVar, zeta, activation and xi. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static Constraint of the SatelliteBlock
 /** Generates the constraints (1), (2) and (3) of the class comments and
  * the bounds on the threshold. */

 void generate_abstract_constraints(
  Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the Objective of the SatelliteBlock, i.e., min zeta

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*----------- Methods for reading the data of the SatelliteBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the SatelliteBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of revisit periods of target k

 Index get_period( Index k ) const { return( periods[ k ] ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of candidate orbits

 Index get_numOrbit( void ) const { return( OrbitSet ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of targets

 Index get_n( void ) const { return( n ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of time stamps

 Index get_t( void ) const { return( t ); }

/*--------------------------------------------------------------------------*/
 /// returns the time horizon [s]

 double get_horizon( void ) const { return( T ); }

/*--------------------------------------------------------------------------*/
 /// returns the time step [s]

 double get_timeStep( void ) const { return( dt ); }

/*--------------------------------------------------------------------------*/
 /// returns the maximum threshold \f$ \theta^{\max} \f$

 double get_theta( void ) const { return( thetaVal ); }

/*--------------------------------------------------------------------------*/
 /// returns the half-aperture of the sensor

 double get_alpha( void ) const { return( alphaHalf ); }

/*--------------------------------------------------------------------------*/
 /// returns the altitude of the orbits [m]

 double get_alt( void ) const { return( altitudeVal ); }

/*--------------------------------------------------------------------------*/
 /// returns the latitude distance of orbit i from target n at time stamp t

 double get_Delta_lat( Index i , Index n , Index t ) const {
  return( CoverageSatLat[ n ][ t ][ i ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the longitude distance of orbit i from target n at time stamp t

 double get_Delta_long( Index i , Index n , Index t ) const {
  return( CoverageSatLong[ n ][ t ][ i ] );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading and writing the Variable --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading and writing the Variable
 *
 * These all require that generate_abstract_variables() has been called.
 *  @{ */

 /// returns a pointer to the Variable xi[ n ][ t ]

 [[nodiscard]] ColVariable * i2p_r( Index n , Index t ) const {
  return( const_cast< ColVariable * >( &xi[ n ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the Variable zeta

 [[nodiscard]] ColVariable * i2p_z( void ) const {
  return( const_cast< ColVariable * >( &zeta[ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the Variable thetaVar

 [[nodiscard]] ColVariable * i2p_theta( void ) const {
  return( const_cast< ColVariable * >( &thetaVar[ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable zeta

 FNumber get_zeta( void ) const { return( zeta[ 0 ].get_value() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable thetaVar

 double get_thetaVar( void ) const { return( thetaVar[ 0 ].get_value() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable xi[ n ][ t ]

 FNumber get_xi( Index n , Index t ) const {
  return( xi[ n ][ t ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current value of the Variable activation[ j ]

 FNumber get_activation( Index j ) const {
  return( activation[ j ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// sets the value of zeta out of the range rng of the values in fstrt
 /** Sets the value of the Variable zeta to *fstrt if rng.first == 0 (zeta
  * being a single Variable, there is nothing to do otherwise). */

 void set_zeta( c_Vec_FNumber_it fstrt ,
                Range rng = Range( 0 , Inf< Index >() ) );

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable xi[ n ][ t ]

 void set_xi1( Index n , Index t , int value ) {
  xi[ n ][ t ].set_value( value );
  }

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable activation[ j ]

 void set_activation1( Index j , int value ) {
  activation[ j ].set_value( value );
  }

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable zeta

 void set_zeta1( int value ) { zeta[ 0 ].set_value( value ); }

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable thetaVar

 void set_thetaVar1( double value ) { thetaVar[ 0 ].set_value( value ); }

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable activation[ j ] and fixes it

 void set_activation( Index j , int value ) {
  activation[ j ].set_value( value );
  activation[ j ].is_fixed( true );
  }

/*--------------------------------------------------------------------------*/
 /// sets the value of the Variable zeta and fixes it

 void set_zeta( int value ) {
  zeta[ 0 ].set_value( value );
  zeta[ 0 ].is_fixed( true );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the
  * Variable of the SatelliteBlock is approximately feasible, i.e., every
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
 /// optimality of a SatelliteBlock is not checked: returns false
 /** Deciding whether the current solution is optimal would require to
  * solve the SatelliteBlock, which is a MILP: this is left to the Solver
  * [see SatelliteSolver.h], and this method always returns false. */

 bool is_optimal( bool useabstract = false ,
                  Configuration * optc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a SatelliteSolution with the current value of zeta
 /** Returns a SatelliteSolution sized for the Variable zeta of this
  * SatelliteBlock; unless emptys == true, it also holds its current value.
  * The Configuration solc is ignored. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the SatelliteBlock
 /** Method for handling Modification. Any "abstract Modification" (one for
  * which Modification::concerns_Block() is true) has its concerns_Block()
  * value set to false and is passed to guts_of_add_Modification() before
  * being forwarded to Block::add_Modification(). Since the data of the
  * SatelliteBlock are not changed by its abstract representation, the
  * abstract Modification require no translation and
  * guts_of_add_Modification() does nothing; changes of the coefficients of
  * the Objective (e.g., those of a Lagrangian relaxation) are just
  * forwarded to the Solver. */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR PRINTING & SAVING THE SatelliteBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the SatelliteBlock
 *  @{ */

 /// print the SatelliteBlock on an ostream
 /** Prints the size of the SatelliteBlock (targets, time stamps, candidate
  * orbits) and its maximum threshold; vlvl is ignored. */

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

 Index n;        ///< the number of targets
 Index t;        ///< the number of time stamps
 Index OrbitSet; ///< the number of candidate orbits

 FNumber dt;          ///< the time step [s]
 FNumber T;           ///< the time horizon [s]
 FNumber alphaHalf;   ///< the half-aperture of the sensor
 FNumber altitudeVal; ///< the altitude of the orbits [m]
 FNumber thetaVal;    ///< the maximum threshold theta^{\max}

 Vec_CNumber periods; ///< the number of revisit periods per target

 boost::multi_array< double , 3 > CoverageSatLat;
 ///< CoverageSatLat[ m ][ t ][ c ], the latitude distances
 boost::multi_array< double , 3 > CoverageSatLong;
 ///< CoverageSatLong[ m ][ t ][ c ], the longitude distances

 std::vector< ColVariable > thetaVar;   ///< the threshold (one Variable)
 std::vector< ColVariable > zeta;       ///< the activation (one Variable)
 std::vector< ColVariable > activation; ///< the orbit Variable pi[ c ]
 boost::multi_array< ColVariable , 2 > xi;
 ///< the observation Variable xi[ m ][ t ]

 std::vector< FRowConstraint > orbitSelection; ///< the constraint (2)
 std::vector< FRowConstraint > thetaUB;        ///< theta <= theta^{\max} zeta
 std::vector< FRowConstraint > thetaLB; ///< theta >= theta^{\max} zeta / 20
 boost::multi_array< FRowConstraint , 2 > activationSat_cnst;
 ///< the constraints xi[ m ][ t ] <= zeta (3)
 std::vector< FRowConstraint > activationSat_cnst_1;
 ///< the constraint zeta <= sum_{ m , t } xi (3)
 boost::multi_array< FRowConstraint , 2 > obs2_cnst;
 ///< the constraints (1) for the latitude
 boost::multi_array< FRowConstraint , 2 > obs4_cnst;
 ///< the constraints (1) for the longitude

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

 SMSpp_insert_in_factory_h; // insert SatelliteBlock in the Block factory

/*--------------------------------------------------------------------------*/

 }; // end( class( SatelliteBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SatelliteSolution -------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a SatelliteBlock
/** A SatelliteSolution holds the value of the Variable zeta of a
 * SatelliteBlock, i.e., whether the satellite is active. An empty
 * SatelliteSolution (one with no value) reads and writes nothing. */

class SatelliteSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*-------------------------------- FRIENDS ---------------------------------*/

 friend SatelliteBlock; ///< make SatelliteBlock friend

/*------------- CONSTRUCTING AND DESTRUCTING SatelliteSolution -------------*/

 /// constructor, it has nothing to do

 explicit SatelliteSolution( void ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// de-serialize a SatelliteSolution out of netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// destructor: it is virtual, and empty

 ~SatelliteSolution() override = default;

/*--------- METHODS DESCRIBING THE BEHAVIOR OF A SatelliteSolution ---------*/

 /// read the value of zeta out of the given SatelliteBlock
 /** Reads the current value of the Variable zeta of the given
  * SatelliteBlock; does nothing if the SatelliteSolution is empty. */

 void read( const Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// write the value of zeta into the given SatelliteBlock
 /** Writes the stored value into the Variable zeta of the given
  * SatelliteBlock; does nothing if the SatelliteSolution is empty. */

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// serialize a SatelliteSolution into a netCDF::NcGroup
 /** Not implemented: exception is thrown. */

 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a SatelliteSolution with the value multiplied by factor

 SatelliteSolution * scale( double factor ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// adds the given SatelliteSolution, multiplied by multiplier
 /** Adds multiplier times the value of the given SatelliteSolution, which
  * must have the same size, to the stored one. */

 void sum( const Solution * solution , double multiplier ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns a copy of the SatelliteSolution
 /** Returns a copy of the SatelliteSolution; if empty == true the copy has
  * the same size but its value is zero. */

 SatelliteSolution * clone( bool empty = false ) const override final;

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// print the SatelliteSolution

 void print( std::ostream & output ) const override final;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*----------------------------- PRIVATE FIELDS -----------------------------*/

 SatelliteBlock::Vec_FNumber v_zeta; ///< the value of zeta, if any

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( SatelliteSolution ) )

/** @} end( group( SatelliteBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SatelliteBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File SatelliteBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
