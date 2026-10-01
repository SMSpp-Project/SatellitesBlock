/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteSatelliteBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DiscreteSatelliteBlock and
 * DiscreteSatelliteSolution classes.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteSatelliteBlock.h"

#include <algorithm>

#include <cmath>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

using FNumber = DiscreteSatelliteBlock::FNumber;
using CNumber = DiscreteSatelliteBlock::CNumber;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DiscreteSatelliteBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( DiscreteSatelliteBlock );

// register DiscreteSatelliteSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( DiscreteSatelliteSolution );

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF DiscreteSatelliteBlock --------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::load( Index n_orbits , Index n_levels )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 OrbitSet = n_orbits;
 ell = n_levels;

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( DiscreteSatelliteBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::load( std::istream & input , char frmt )
{
 Index n_orbits;
 Index n_levels;
 if( ! ( input >> eatcomments >> n_orbits ) )
  throw( std::invalid_argument(
   "DiscreteSatelliteBlock::load: error reading the number of orbits" ) );

 if( ! ( input >> eatcomments >> n_levels ) )
  throw( std::invalid_argument(
   "DiscreteSatelliteBlock::load: error reading the number of levels" ) );

 load( n_orbits , n_levels );

 } // end( DiscreteSatelliteBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::generate_abstract_variables(
 Configuration * stvv )
{
 if( AR3 & HasVar ) // the variables are there already
  return;           // nothing to do

 y.resize( boost::extents[ OrbitSet ][ ell ] );
 for( Index o = 0 ; o < OrbitSet ; ++o )
  for( Index l = 0 ; l < ell ; ++l )
   y[ o ][ l ].set_type( ColVariable::kBinary , eNoMod );

 add_static_variable( y , "y" );

 AR3 |= HasVar;

 } // end( DiscreteSatelliteBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::generate_abstract_constraints(
 Configuration * stcc )
{
 if( AR2 & HasCnst ) // the constraints are there already
  return;            // nothing to do

 LinearFunction::v_coeff_pair v_var;
 v_var.reserve( OrbitSet * ell );
 for( Index o = 0 ; o < OrbitSet ; ++o )
  for( Index l = 0 ; l < ell ; ++l )
   v_var.push_back( std::make_pair( &y[ o ][ l ] , 1.0 ) );

 one.set_function( new LinearFunction( std::move( v_var ) ) , eNoMod );
 one.set_rhs( 1.0 , eNoMod );
 one.set_lhs( -Inf< double >() , eNoMod );

 add_static_constraint( one , "one" );

 AR2 |= HasCnst;

 } // end( DiscreteSatelliteBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::generate_objective( Configuration * objc )
{
 if( AR1 & HasObj ) // the objective is there already
  return;           // cowardly (and silently) return

 LinearFunction::v_coeff_pair v_var;
 v_var.reserve( OrbitSet * ell );
 for( Index o = 0 ; o < OrbitSet ; ++o )
  for( Index l = 0 ; l < ell ; ++l )
   v_var.push_back( std::make_pair( &y[ o ][ l ] , 1.0 ) );

 c.set_function( new LinearFunction( std::move( v_var ) , 0 ) , eNoMod );
 set_objective( &c , eNoMod );

 AR1 |= HasObj;

 } // end( DiscreteSatelliteBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR CHECKING THE Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool DiscreteSatelliteBlock::is_feasible( bool useabstract ,
                                          Configuration * fsbc )
{
 if( ! ( AR3 & HasVar ) ) // there is no solution to check
  return( false );

 FNumber eps = 0;
 auto tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
   f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 FNumber sum = 0;
 for( Index o = 0 ; o < OrbitSet ; ++o )
  for( Index l = 0 ; l < ell ; ++l ) {
   const FNumber val = y[ o ][ l ].get_value();
   if( ( val < -eps ) || ( val > 1 + eps ) ||
       ( std::abs( val - std::round( val ) ) > eps ) )
    return( false );
   sum += val;
   }

 return( sum <= 1 + eps );

 } // end( DiscreteSatelliteBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool DiscreteSatelliteBlock::is_optimal( bool useabstract ,
                                         Configuration * optc )
{
 if( ! ( ( AR3 & HasVar ) && ( AR1 & HasObj ) ) )
  return( false );

 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  using CFpair = std::pair< CNumber , FNumber >;
  if( auto toptc = dynamic_cast< SimpleConfiguration< CFpair > * >( optc ) ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
   }
  else {
   if( auto ttoptc =
        dynamic_cast< SimpleConfiguration< CNumber > * >( optc ) )
    ceps = ttoptc->f_value;

   if( f_BlockConfig && f_BlockConfig->f_is_feasible_Configuration )
    if( auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
         f_BlockConfig->f_is_feasible_Configuration ) )
     feps = fsbc->f_value;
   }
  }
 else if( f_BlockConfig ) {
  if( f_BlockConfig->f_is_optimal_Configuration )
   if( auto csbc = dynamic_cast< SimpleConfiguration< CNumber > * >(
        f_BlockConfig->f_is_optimal_Configuration ) )
    ceps = csbc->f_value;

  if( f_BlockConfig->f_is_feasible_Configuration )
   if( auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
        f_BlockConfig->f_is_feasible_Configuration ) )
    feps = fsbc->f_value;
  }

 SimpleConfiguration< FNumber > feascfg( feps );
 if( ! is_feasible( useabstract , &feascfg ) )
  return( false );

 auto lf = dynamic_cast< LinearFunction * >( c.get_function() );
 if( ! lf )
  return( false );

 // since at most one y is nonzero, the optimal value is the minimum
 // between 0 and the smallest coefficient (plus the constant term, which
 // is left out of both values)
 CNumber best = 0;
 CNumber value = 0;
 for( const auto & [ var , coeff ] : lf->get_v_var() ) {
  best = std::min( best , coeff );
  value += coeff * var->get_value();
  }

 return( value - best <= ceps * std::max( CNumber( 1 ) , std::abs( best ) ) );

 } // end( DiscreteSatelliteBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/

Solution * DiscreteSatelliteBlock::get_Solution( Configuration * solc ,
                                                 bool emptys )
{
 auto * sol = new DiscreteSatelliteSolution();

 if( AR3 & HasVar ) {
  sol->v_y.assign( OrbitSet * ell , 0 );
  if( ! emptys )
   sol->read( this );
  }

 return( sol );

 } // end( DiscreteSatelliteBlock::get_Solution )

/*--------------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 } // end( DiscreteSatelliteBlock::add_Modification )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE DiscreteSatelliteBlock */
/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::print( std::ostream & output , char vlvl ) const
{
 output << OrbitSet << " " << ell << std::endl;

 } // end( DiscreteSatelliteBlock::print )

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteSatelliteBlock::guts_of_destructor( void )
{
 // clear() the Constraint and the Objective, so that they do not bother to
 // un-register themselves from Variable that are going to be deleted anyway
 one.clear();
 c.clear();

 // explicitly reset all Constraint and Variable, so that a new abstract
 // representation is not added to the (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 y.resize( boost::extents[ 0 ][ 0 ] );

 AR1 = AR2 = AR3 = 0;

 } // end( DiscreteSatelliteBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/
// the DiscreteSatelliteBlock has no physical representation that an
// abstract Modification may change, hence there is nothing to do here

void DiscreteSatelliteBlock::guts_of_add_Modification( p_Mod mod ,
                                                       ChnlName chnl )
{} // end( DiscreteSatelliteBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF DiscreteSatelliteSolution ------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteSatelliteSolution::deserialize( const netCDF::NcGroup & group )
{
 throw( std::logic_error(
  "DiscreteSatelliteSolution::deserialize: not implemented" ) );
 }

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const DiscreteSatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "DiscreteSatelliteSolution::read: block is "
                                "not a DiscreteSatelliteBlock" ) );

 if( v_y.empty() )
  return;

 if( v_y.size() != SATB->y.num_elements() )
  throw( std::invalid_argument(
   "DiscreteSatelliteSolution::read: wrong Block size" ) );

 auto var = SATB->y.data();
 for( auto & val : v_y )
  val = ( var++ )->get_value();

 } // end( DiscreteSatelliteSolution::read )

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteSolution::write( Block * block )
{
 auto SATB = dynamic_cast< DiscreteSatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "DiscreteSatelliteSolution::write: block is "
                                "not a DiscreteSatelliteBlock" ) );

 if( v_y.empty() )
  return;

 if( v_y.size() != SATB->y.num_elements() )
  throw( std::invalid_argument(
   "DiscreteSatelliteSolution::write: wrong Block size" ) );

 auto var = SATB->y.data();
 for( auto val : v_y )
  ( var++ )->set_value( val );

 } // end( DiscreteSatelliteSolution::write )

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteSolution::serialize( netCDF::NcGroup & group ) const
{
 throw( std::logic_error(
  "DiscreteSatelliteSolution::serialize: not implemented" ) );
 }

/*--------------------------------------------------------------------------*/

DiscreteSatelliteSolution * DiscreteSatelliteSolution::scale(
 double factor ) const
{
 auto * sol = clone();
 for( auto & val : sol->v_y )
  val *= factor;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteSolution::sum( const Solution * solution ,
                                     double multiplier )
{
 auto tsol = dynamic_cast< const DiscreteSatelliteSolution * >( solution );
 if( ! tsol )
  throw( std::invalid_argument( "DiscreteSatelliteSolution::sum: solution is "
                                "not a DiscreteSatelliteSolution" ) );

 if( tsol->v_y.size() != v_y.size() )
  throw(
   std::invalid_argument( "DiscreteSatelliteSolution::sum: wrong size" ) );

 for( Index i = 0 ; i < v_y.size() ; ++i )
  v_y[ i ] += multiplier * tsol->v_y[ i ];
 }

/*--------------------------------------------------------------------------*/

DiscreteSatelliteSolution * DiscreteSatelliteSolution::clone(
 bool empty ) const
{
 auto * sol = new DiscreteSatelliteSolution();

 if( empty )
  sol->v_y.assign( v_y.size() , 0 );
 else
  sol->v_y = v_y;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void DiscreteSatelliteSolution::print( std::ostream & output ) const
{
 output << "DiscreteSatelliteSolution [" << v_y.size() << "]:";
 for( auto val : v_y )
  output << " " << val;
 output << std::endl;
 }

/*--------------------------------------------------------------------------*/
/*------------------ End File DiscreteSatelliteBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
