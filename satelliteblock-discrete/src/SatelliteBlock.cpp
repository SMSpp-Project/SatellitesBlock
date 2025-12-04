/*--------------------------------------------------------------------------*/
/*--------------------- File Satellite.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SatelliteBlock class.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SatelliteBlock.h"
#include "DQuadFunction.h"
#include <iomanip>
#include <cmath>

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures representing the
  * abstract and the physical representations agree. */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using FNumber = SatelliteBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const double RAYON = 6378136.3;
static const double M_limit = 10*3.14159265;
static const auto PI = 3.14159265;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

static void print_UB( std::ostream & os , FNumber ub )
{
 if( ub == Inf< FNumber >() )
  os << "+Inf";
 else
  os << ub;
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ

template< typename T >
static Index countdiff( T beg , T end , T cmp )
{
 Index ndiff = 0;
 for( ; beg != end ; )
  if( *(beg++) != *(cmp++) )
   ndiff++;

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector< T > & vec , c_Subset & nms ,
		      typename std::vector< T >::const_iterator cmp ,
		      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices

template< typename T >
static Index countdiff( std::vector< T > & vec , c_Subset & nms ,
			typename std::vector< T >::const_iterator cmp ,
			Index n_max )
{
 Index ndiff = 0;
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   ndiff++;
  }

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another

template< typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SatelliteBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SatelliteBlock );

// register SatelliteSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( SatelliteSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF SatelliteBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( FNumber indOrbit , FNumber indTheta )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

  ell = indTheta;
  OrbitSet = indOrbit;

 // allocate observability variables - - - - - - - - - - - - - - - - - - - - - 

 //generate_abstract_variables();
 //generate_objective();

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 //if( anyone_there() )
 // add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( SatelliteBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 // allocate observability variables - - - - - - - - - - - - - - - - - - - - - 

 //generate_abstract_variables();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( SatelliteBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_abstract_variables( Configuration *stvv )
{

  if( AR3 & HasVar )  // the variables are there already
    return;           // nothing to do

  y.resize( boost::extents[ OrbitSet ][ ell ] );
  for( Index i = 0 ; i < OrbitSet ; ++i )
    for( Index j = 0 ; j < ell ; ++j )
      y[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( y );

  AR3 |= HasVar;

 }  // end( SatelliteBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_abstract_constraints( Configuration *stcc )
{

  if( AR2 & HasCnst )  // the constraints are there already
    return;           // nothing to do
  
  LinearFunction::v_coeff_pair v_var;
  for( Index i = 0 ; i < OrbitSet ; ++i ){
    for( Index j = 0 ; j < ell ; ++j ){
      v_var.push_back( std::make_pair( &y[ i ][ j ], 1.0 ));
    }
  }

  LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
  one.set_rhs( 1.0 );
  one.set_lhs( -Inf< double >() ); 
  one.set_function( FunctSat );
  add_static_constraint( one );

  AR2 |= HasCnst;

 }  // end( SatelliteBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_objective( Configuration *objc )
{

 if( AR1 & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

  DQuadFunction::v_coeff_triple p;
  LinearFunction::v_coeff_pair P;

  for( Index i = 0 ; i < OrbitSet ; ++i ){
    for( Index j = 0 ; j < ell ; ++j ){
        p.push_back( std::make_tuple( &y[ i ][ j ], 1.0, 0.0 ));
        P.push_back( std::make_pair( &y[ i ][ j ], 1.0 ));
    }
  }
  
  //c.set_function( new DQuadFunction( std::move( p ) , 0 ) , eNoMod );
  c.set_function( new LinearFunction( std::move( P ) , 0 ) , eNoMod );
  set_objective( & c , eNoMod );

  AR1 |= HasObj;

 }  // end( SatelliteBlock::generate_objective )

/*--------------------------------------------------------------------------*/

 bool SatelliteBlock::is_feasible( bool useabstract , Configuration *fsbc )
{
 FNumber eps = 0;
 auto tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                        f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return( 0 );

 }  // end( SatelliteBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool SatelliteBlock::is_optimal( bool useabstract , Configuration *optc )
{
 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  if( auto toptc =
      dynamic_cast< SimpleConfiguration< std::pair< CNumber , FNumber > > * >(
								    optc ) ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
   }
  else {
   auto ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >( optc );

   if( ( ! ttoptc ) && f_BlockConfig &&
       f_BlockConfig->f_is_optimal_Configuration )
    ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                           f_BlockConfig->f_is_optimal_Configuration );
   if( ttoptc )
    ceps = ttoptc->f_value;

   if( f_BlockConfig && f_BlockConfig->f_is_feasible_Configuration ) {
    auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration );
    if( fsbc )
     feps = fsbc->f_value;
    }
   }
  }
 else
  if( f_BlockConfig ) {
   if( f_BlockConfig->f_is_optimal_Configuration )
    if( auto csbc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                              f_BlockConfig->f_is_optimal_Configuration ) )
     ceps = csbc->f_value;

   if( f_BlockConfig->f_is_feasible_Configuration )
    if( auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration ) )
     feps = fsbc->f_value;
   }

 return( 0 );

 }  //  end( SatelliteBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/


 Solution * SatelliteBlock::get_Solution( Configuration * solc , bool emptys )
{
   
 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;
 
 auto *sol = new SatelliteSolution();

 if( ! emptys )
  sol->read( this );
 
 return( sol );

 }  // end( SatelliteBlock::get_Solution )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE SatelliteBlock ---*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::print( std::ostream  & output , char vlvl ) const
{
 
 }  // end( SatelliteBlock::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints

 c.clear();  // clear the Objective

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 //reset_dynamic_constraints();
 //reset_dynamic_variables();
 reset_objective();

 }  // end( SatelliteBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since SatelliteBlock is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN SatelliteBlock WHEN THIS METHOD IS
  *   EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS ISSUED:
  *   NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some of the logic here.*/

  // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {

  auto lfo = static_cast< LinearFunction * const >( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Range in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble

  return;
  }
*/

 //throw( std::invalid_argument( "unsupported Modification to SatelliteBlock" ) );

 }  // end( SatelliteBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF DCRSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteSolution::deserialize( const netCDF::NcGroup & group )
{}

void SatelliteSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const SatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a SatelliteBlock" ) );

}

void SatelliteSolution::write( Block * block ) 
{

 auto SATB = dynamic_cast<SatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a SatelliteBlock" ) );

}

void SatelliteSolution::serialize( netCDF::NcGroup & group ) const
{}

SatelliteSolution * SatelliteSolution::scale( double factor ) const
{
  auto * sol = SatelliteSolution::clone( true );
  return( sol );
}

void SatelliteSolution::sum( const Solution * solution , double multiplier )
{}

SatelliteSolution * SatelliteSolution::clone( bool empty ) const
{
  auto * sol = new SatelliteSolution();
 
 return( sol );
}

/*--------------------------------------------------------------------------*/
/*------------------- End File SatelliteBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
