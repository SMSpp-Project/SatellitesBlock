/*--------------------------------------------------------------------------*/
/*------------------------- File DiscreteConstellationBlock.cpp ------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DiscreteConstellationBlock class.
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

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteConstellationBlock.h"

#include <math.h>

#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC MEMBERS -------------------------------*/
/*--------------------------------------------------------------------------*/

// register DiscreteConstellationBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( DiscreteConstellationBlock );

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const auto RAYON = 6378136.3; //[m]
static const auto PI = 3.14159265;
static const auto MU = 3.986004418e14; //[m^3/s^2]
static const auto WE = 7.2921e-5;
static const auto FACTOR = 1.2;
static const auto angle0 = -1.3882860164509252;

/*--------------------------------------------------------------------------*/
/*------------------------ OTHER INITIALIZATIONS ---------------------------*/
/*--------------------------------------------------------------------------*/

void  DiscreteConstellationBlock::load( std::istream & input , char frmt ){
}

/*--------------------------------------------------------------------------*/
// this load() is the discretized-theta counterpart of
// ConstellationBlock::load( const std::string & ) [see there for a full
// account of the orbital-mechanics computations, which are identical here
// up to and including the discretization of inclination/RAAN/mean anomaly
// and the ground-track propagation]: the only substantial addition is that,
// besides the candidate orbits [C], the observability threshold theta^{\max}
// is here also discretized into ell (== 3) decreasing levels thetaValF[ l ],
// and, for every satellite / time stamp / target / orbit / level quintuple,
// the corresponding 0/1 observability outcome is precomputed once and for
// all into obs[][][][][], to be used directly (as a constraint coefficient,
// not as a big-M linearization) by generate_abstract_constraints() below

void DiscreteConstellationBlock::load( const std::string & input , char frmt )
{

 // ensure starting from clean slate
 guts_of_destructor();

 std::ifstream iFile( input );

 iFile >> horizon;
 horizon *= 3600.0;
 iFile >> time_step;
 iFile >> targets;
  

 std::cout << "TARGETS: " << targets << "\n";
 std::cout << "TIME STEP: " << time_step << "\n";

 Vec_CNumber Latitude;
 Vec_CNumber Longitude;

 double altitudeSetVal;
 Vec_CNumber periodSat;
 int indexLen = 0;
 double j1 = 0.0;

 for( Index i = 0 ; i < horizon/3600.0 ; ++i ){
   j1++;
	altitudeSetVal = cbrt( ( MU* pow((horizon)/j1,2.0) ) / (4.0*pow(PI,2.0))) - RAYON;
	if ( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 )
   {
		indexLen++;
   }
 }

 FNumber altSet = indexLen;
 altitude.resize( indexLen );
 periodSat.resize( indexLen );
 indexLen = 0;

 j1 = 0.0;
 for( Index i = 0 ; i < horizon/3600.0 ; ++i )
 {
   j1++;
	altitudeSetVal = cbrt( ( MU * pow((horizon)/j1, 2.0) ) / (4.0*pow(PI, 2.0) )) - RAYON;
   //std::cout << altitudeSetVal << "\n";
	if ( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 )
   {
		altitude[indexLen] = altitudeSetVal;
		periodSat[indexLen] = horizon/j1;
      indexLen++;
   }
 }

 Latitude.resize( targets );
 Longitude.resize( targets );
 periods.resize( targets );
 for( Index i = 0 ; i < targets ; ++i )
 {
   iFile >> periods[i];
   iFile >> Latitude[i];
   iFile >> Longitude[i];
   Latitude[i] *= PI/180;
   Longitude[i] *= PI/180;
 }

 Vec_CNumber t_p;
 Vec_CNumber t_u;
 Vec_CNumber t_GM;
 Vec_CNumber thetaVal;

 t_p.resize(altSet);
 t_u.resize(altSet);
 t_GM.resize(altSet);
 thetaVal.resize(altSet);

 double alphalim;
 double Theta_min = ((2*PI*time_step)/(2*periodSat[indexLen-1]))*FACTOR;
 std::cout << "theta_min: " << Theta_min << "\n";
 double aHalf = atan(sin(Theta_min)/((RAYON+altitude[indexLen-1])/RAYON - cos(Theta_min)));

 for( Index ii = 0 ; ii < altSet ; ++ii ){
   if(((RAYON+altitude[ ii ])/RAYON)*sin(aHalf) > 1){
      aHalf = asin((RAYON/(RAYON+altitude[ ii ])));
      std::cout << "WARNING: computed alpha_lim\n";
      break;
   } 
 }

 //aHalf = 50 * PI/180;

 double alt = altitude[altSet-1];
 ///double alt = altitude[0];
 for( Index ii = 0 ; ii < altSet ; ++ii ){
   t_p[ ii ] = sqrt(MU / ( pow(RAYON+altitude[ii],3)));
   t_u[ ii ] = sqrt(MU / (RAYON+altitude[ ii ]));
   t_GM[ ii ] = sqrt((RAYON+altitude[ ii ]) / MU);
   thetaVal[ ii ] = -aHalf + asin(((RAYON+altitude[ii])/RAYON)*sin(aHalf));
 }

 double Theta_max = thetaVal[altSet-1];
 ///Theta_min = thetaVal[0];
 //Theta_min = thetaVal[altSet-1];
 std::cout << "theta_min: " << Theta_min << "\n";
 Theta_max = Theta_min;
 double numbOfDiscretize = ceil(PI/Theta_min);
 
 std::cout << "numbOfDiscretize: " << numbOfDiscretize << "\n";

 FNumber incSet = numbOfDiscretize;
 FNumber ascSet = numbOfDiscretize;
 FNumber anmSet = numbOfDiscretize;

 double start_in = 0;
 double end_in = PI;
 size_t num_in = numbOfDiscretize;

 double dx = (end_in - start_in) / (num_in - 1);
 std::vector<double> x(num_in);
 int iter = 0;
 std::generate(x.begin(), x.end(), [&] { return start_in + (iter++) * dx; });

 Vec_CNumber inclination = x;
 //std::cout << inclination;

 end_in = 2*PI;
 dx = (end_in - start_in) / (num_in - 1);
 iter = 0;
 std::generate(x.begin(), x.end(), [&] { return start_in + (iter++) * dx; });

 Vec_CNumber nodeAscendant = x;
 Vec_CNumber meanAnomaly = x;

 double thetaValFinal = Theta_min;
 double altitudeFinal = altitude[indexLen-1]; 
 ///double altitudeFinal = altitude[0];
 
 //std::cout << "theta_min: " << Theta_min << "\n";
 std::cout << "alpha_half: " << aHalf*180/PI << "\n";
 std::cout << "theta: " << thetaValFinal << "\n";
 std::cout << "altitude: " << altitudeFinal << "\n";

 FNumber t = horizon / time_step;

 boost::multi_array< double , 3 > CoverageSatLat(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLat1(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong1(boost::extents[targets][t][incSet*ascSet*anmSet]);

 // ell discretized observability-threshold levels are used for every
 // satellite; thetaValF[] is (re)computed per-satellite a few lines below

 int ell = 3;
 thetaValF.resize(ell);

 double lat_Sat;
 double long_Sat;

 int indexOrbit;
 int indexOrbit1;

 iFile >> satellites;
 obs.resize(boost::extents[satellites][t][targets][incSet*ascSet*anmSet][ell]);
 //std::cout << "number Orbits: " << indexOrbit << "\n";
 /*
 satellites = 0;
 for( Index i = 0 ; i < targets ; ++i )
 {
   satellites += periods[i];
 }
 */
 v_Block.resize( satellites );
 indexOrbitSat.resize( satellites );
 ThetaValMaxSat.resize( satellites );
 thetaValSat.resize( boost::extents[satellites][ell] );
 ellSat.resize( satellites );
 std::cout << "SATELLITES: " << satellites << "\n";

for(Index isat = 0; isat < satellites; ++isat){
 //int ii = indexLen-1;

 int ialt = std::floor(satellites/3);
 int ii = 0;
 if( isat < ialt )
   ii = 2;
 if( isat > 2*ialt)
   ii = 0;
 if(isat >= ialt and isat <= 2*ialt)
   ii = 1;

 std::cout << altitude[ii] << std::endl;

 ///int ii = 0; 
 int index1 = -1;
 indexOrbit = 0;
 int addOrbit = 0;
 indexOrbit1 = 0;
 //double thetaValF = FACTOR*(satellites-isat)/satellites * thetaValFinal;
 //double thetaValF = (satellites - isat)/satellites * thetaValFinal;
 // level 0 is the loosest (largest) threshold, level ell-1 the tightest;
 // note that only thetaValF[ 0 ] is actually used below to decide whether
 // an orbit is worth keeping at all (indexOrbit1), while all ell levels
 // are used afterwards to fill in obs[][][][][]

 for(double l = 0; l < ell; l++){
   thetaValF[l] = (ell-l)/ell * thetaValFinal/10.0;
   std::cout << "theta=" << thetaValF[l] << std::endl;
 }

    for( Index jj = 0 ; jj < incSet ; ++jj ) 
    {
       for( Index k = 0 ; k < ascSet ; ++k ) 
       {
         for( Index l = 0 ; l < anmSet ; ++l )
          {
            index1 += 1;
            indexOrbit1 = 0;
            for( Index j = 0 ; j < t ; ++j ) 
            {
             lat_Sat = asin(((sin(inclination[jj])*(altitude[ii]+RAYON)*sin(meanAnomaly[l]))*cos(t_p[ii]*((j)*time_step))/(altitude[ii]+RAYON)) 
		                        + ((sin(inclination[jj])*t_u[ii]*cos(meanAnomaly[l]))*sin(t_p[ii]*((j)*time_step))*t_GM[ii]));
             
             long_Sat = fmod(-(angle0 + (WE*((j)*time_step))) + 
		                        atan2((((sin(nodeAscendant[k])*(altitude[ii]+RAYON)*cos(meanAnomaly[l])) + (cos(nodeAscendant[k])*cos(inclination[jj])*
		                        (altitude[ii]+RAYON)*sin(meanAnomaly[l])))*cos(t_p[ii]*((j)*time_step))/(altitude[ii]+RAYON))	+ ((-(sin(nodeAscendant[k])*t_u[ii]*sin(meanAnomaly[l]))
		                        + (cos(nodeAscendant[k])*cos(inclination[jj])*t_u[ii]*cos(meanAnomaly[l])))*sin(t_p[ii]*((j)*time_step))*t_GM[ii]),
		                        ((((cos(nodeAscendant[k])*(altitude[ii]+RAYON)*cos(meanAnomaly[l])) - (sin(nodeAscendant[k])*cos(inclination[jj])*
		                        (altitude[ii]+RAYON)*sin(meanAnomaly[l])))* cos(t_p[ii]*(j*time_step))/(altitude[ii]+RAYON)) + ((-(cos(nodeAscendant[k])*t_u[ii]*sin(meanAnomaly[l]))
		                        -(sin(nodeAscendant[k])*cos(inclination[jj])*t_u[ii]*cos(meanAnomaly[l])))*sin((t_p[ii]*(j*time_step)))*t_GM[ii]))), (2*PI));
             if( long_Sat <= 0 )
			      long_Sat += 2*PI;

	     /*
	     for( Index i = 0 ; i < targets ; ++i )                                                                                                                                                        
             {                                                                                                                                                                                              
                CoverageSatLat[ i ][ j ][ index1 ] = std::abs(Latitude[i] - lat_Sat);                                                                                                                      
                CoverageSatLong[ i ][ j ][ index1 ] = std::abs(Longitude[i] - long_Sat) * cos(Latitude[i]);
	     }}
	     indexOrbit += 1;
	     */
	     for( Index i = 0 ; i < targets ; ++i ) 
             {
                CoverageSatLat1[ i ][ j ][ index1 ] = std::abs(2 * asin( 0.5 * sqrt(1 - cos(Latitude[i]-lat_Sat))));
                //std::abs(Latitude[i] - lat_Sat); 
                CoverageSatLong1[ i ][ j ][ index1 ] = std::abs(2 * asin( 0.5 * sqrt(cos(Latitude[i]) * cos(Latitude[i]) * (1 - cos(Longitude[i] - long_Sat))))) * cos(Latitude[i]);
                //std::abs(Longitude[i] - long_Sat) * cos(Latitude[i]); 
                if(cos(Latitude[i]) < 0)
                  std::cout << "ERROR!" << "\n";
               //std::cout << CoverageSatLat1[ i ][ j ][ index1 ] << " " << CoverageSatLong1[ i ][ j ][ index1 ] << std::endl;
                //std::cout << i << " " << j << " " << ii << " " << jj << " " << k << " " << l << " " << "\n";
                if (CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValF[0] and CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValF[0])
                  indexOrbit1 += 1;
             }
            }
            // as in ConstellationBlock::load(), an orbit is only kept if it
            // observes at least one target at some time stamp (indexOrbit1
            // uses the loosest threshold thetaValF[ 0 ]); for every kept
            // orbit, obs[ isat ][ j ][ i ][ indexOrbit ][ l1 ] is set to 1
            // iff. the coverage distance is within threshold level l1, for
            // every one of the ell levels: this is the precomputed
            // observability outcome used directly, as a coefficient, by
            // generate_abstract_constraints() below

            if (indexOrbit1>=1){
               for( Index j = 0 ; j < t ; ++j )
               {
                  for( Index i = 0 ; i < targets ; ++i )
                  {
                  CoverageSatLat[ i ][ j ][ indexOrbit ] = CoverageSatLat1[ i ][ j ][ index1 ];
                  CoverageSatLong[ i ][ j ][ indexOrbit ] = CoverageSatLong1[ i ][ j ][ index1 ];
                  for(int l1 = 0; l1 < ell; ++l1){
                     if(CoverageSatLong[ i ][ j ][ indexOrbit ] - thetaValF[l1] <= 0.0 and CoverageSatLat[ i ][ j ][ indexOrbit ] - thetaValF[l1] <= 0.0){
                        obs[isat][j][i][indexOrbit][l1] = 1.0;
                     } else {
                        obs[isat][j][i][indexOrbit][l1] = 0.0;
                     }
                  }
                }
               }
               indexOrbit += 1;
            }
          }
       }
    }

      // record, for this satellite, how many orbits survived (indexOrbit)
      // and its threshold levels (thetaValSat[]), then create and load()
      // its DiscreteSatelliteBlock with the corresponding (orbit, level)
      // grid size; note the ground-track coverage itself (obs[][][][][])
      // is not passed to the DiscreteSatelliteBlock, but kept here and
      // used directly when building the observability constraints below

      indexOrbitSat[isat] = indexOrbit;
      for (Index lll = 0; lll < ell; ++lll)
         thetaValSat[isat][lll] = thetaValF[lll];
      ellSat[isat] = ell;
      std::cout << "number Orbits: " << indexOrbit << " " << ell << "\n";
      auto SB = new DiscreteSatelliteBlock( this );
      SB->load( indexOrbit , ell );
      v_Block[ isat ] = SB;
 }
    
 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( DiscreteConstellationBlock::load( const std::string & input )

/*--------------------------------------------------------------------------*/

// simply delegates to each nested DiscreteSatelliteBlock

void DiscreteConstellationBlock::generate_abstract_variables( Configuration * stvv )
{

  for( auto blck : v_Block )
    blck->generate_abstract_variables();
 }

/*--------------------------------------------------------------------------*/
// first delegates to each nested DiscreteSatelliteBlock (its own
// "at most one (orbit,level)" constraint), then adds the
// DiscreteConstellationBlock-level constraints linking them together:
// thetaM, observation / observation1 (as in ConstellationBlock, but
// expressed via y[][] and the precomputed obs[][][][][] coefficients
// rather than via xi[][] and a big-M linearization), and symmetry breaking

void DiscreteConstellationBlock::generate_abstract_constraints( Configuration * stcc )
{

 if( ! ( AR ) ) {
   for( auto blck : v_Block )
    blck->generate_abstract_constraints();

   // thetaM: for every satellite k and every selected (orbit,level) pair
   // y[ o ][ l ], penalizes/forbids levels l whose threshold thetaValSat[k][l]
   // exceeds 0.9 times the loosest level thetaValSat[k][0], mirroring the
   // per-satellite thetaM bound of ConstellationBlock

   LinearFunction::v_coeff_pair v_var12;
   for( Index k = 0 ; k < satellites ; ++k )
      for( Index o = 0 ; o < indexOrbitSat[k] ; ++o )
         for( Index l = 0 ; l < ellSat[k] ; ++l )
         v_var12.push_back( std::make_pair( static_cast< DiscreteSatelliteBlock * >(
                  v_Block[ k ] )->i2p_y(o,l), thetaValSat[k][l] - 0.9 * thetaValSat[k][0] ));

   thetaM.set_function( new LinearFunction( std::move( v_var12 )));
   thetaM.set_rhs( 0.0 );
   thetaM.set_lhs( -Inf< double >() );

   add_static_constraint( thetaM , "thetaM" );

  // generate the observability constraints  - - - - - - - - - - - - - - -
  // same revisit-time semantics as ConstellationBlock::observation, but
  // the coefficient of y[ o ][ l ] for satellite k, target i and revisit
  // window j is obb = sum_{tt in window} obs[ k ][ tt ][ i ][ o ][ l ],
  // i.e., how many time stamps of that window orbit/level (o,l) of
  // satellite k would observe target i

   double maxPeriods = *max_element(periods.begin(), periods.end());
   observation.resize(boost::multi_array< FRowConstraint , 2 >::extent_gen()[ targets ][ maxPeriods ]);
   double pp;
   double obb;

   for( Index i = 0 ; i < targets ; ++i ) {
    pp = horizon/time_step/periods[i];
    //std::cout << maxPeriods << "\n";
    for( Index j = 0 ; j < periods[i] ; ++j ) {
       LinearFunction::v_coeff_pair v_var;

      for( Index k = 0 ; k < satellites ; ++k ) {
         for( Index o = 0 ; o < indexOrbitSat[k] ; ++o ) {
            for( Index l = 0 ; l < ellSat[k] ; ++l ) {
               obb = 0.0;
                  for( Index tt = j*pp ; tt < (j+1)*pp ; ++tt ) {
                     obb += obs[k][tt][i][o][l];
                  }
               v_var.push_back( std::make_pair( static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->i2p_y( o , l ), obb ));
            }
         }
      }

      // fake constraints when j >= periods[i]
      observation[i][j].set_function( new LinearFunction( std::move( v_var )));
      observation[i][j].set_rhs( Inf< double >() ); 
      observation[i][j].set_lhs( 1.0 );
    }

    for( Index j = periods[i] ; j < maxPeriods ; ++j ) {
       LinearFunction::v_coeff_pair v_var;

    for( Index k = 0 ; k < satellites ; ++k ) 
         for( Index o = 0 ; o < indexOrbitSat[k] ; ++o ) 
            for( Index l = 0 ; l < ellSat[k] ; ++l ) 
               for( Index tt = j*pp ; tt < (j+1)*pp ; ++tt ) 
                  v_var.push_back( std::make_pair( static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->i2p_y( o , l ), 0.0 ));

      // fake constraints when j >= periods[i]
      observation[i][j].set_function( new LinearFunction( std::move( v_var )));
      observation[i][j].set_rhs( 0.0 ); 
      observation[i][j].set_lhs( 0.0 );
    }

   }

   add_static_constraint( observation , "observation" );

   // observation1: per-time-stamp version, mirroring ConstellationBlock's
   // observation1 (at most one satellite observes each target at each
   // time stamp), again expressed via the obs[][][][][] coefficients

   observation1.resize(boost::multi_array< FRowConstraint , 2 >::extent_gen()[ targets ][ horizon/time_step ]);

   for( Index i = 0 ; i < targets ; ++i ) {
    for( Index j = 0 ; j < horizon/time_step ; ++j ) {
       LinearFunction::v_coeff_pair v_var;

   for( Index k = 0 ; k < satellites ; ++k ) {
      for( Index o = 0 ; o < indexOrbitSat[k] ; ++o ) {
         for( Index l = 0 ; l < ellSat[k] ; ++l ) {
            v_var.push_back( std::make_pair( static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->i2p_y( o , l ), obs[k][j][i][o][l] ));
         }
      }
   }

   observation1[i][j].set_function( new LinearFunction( std::move( v_var )));
   observation1[i][j].set_rhs( 1.0 ); 
   observation1[i][j].set_lhs( -Inf< double >() );


   }}

   add_static_constraint( observation1 , "observation1" );

   double ellSatMax = *max_element(ellSat.begin(), ellSat.end());
   double indexOrbitSatMax = *max_element(indexOrbitSat.begin(), indexOrbitSat.end());

   int ialt = std::floor(satellites/3);

   // (currently disabled, see add_static_constraint() being commented out
   // below) symmetry-breaking: for every pair of consecutive satellites i,
   // i+1 belonging to the same altitude-band third (cf. the "ii" grouping
   // in load()), forces sum_{o,l} y_i[o][l] <= sum_{o,l} y_{i+1}[o][l], to
   // cut away solutions that only differ by which interchangeable satellite
   // within the group is the active one; note that symmetry[ i ] is left
   // completely unset (no function/rhs/lhs) whenever i and i+1 are *not*
   // in the same group, which would need to be fixed before re-enabling
   // add_static_constraint() below

   symmetry.resize( satellites-1 );

   for( Index i = 0 ; i < satellites-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_vars;
      if(std::floor(i/3)==std::floor((i+1)/3)){
         for( Index o = 0 ; o < indexOrbitSat[i] ; ++o )
            for( Index l = 0 ; l < ellSat[i] ; ++l )
                  v_vars.push_back( std::make_pair( static_cast< DiscreteSatelliteBlock * >( v_Block[ i ] )->i2p_y( o , l ), 1.0 ));
         for( Index o = 0 ; o < indexOrbitSat[i+1] ; ++o )
            for( Index l = 0 ; l < ellSat[i+1] ; ++l )
                  v_vars.push_back( std::make_pair( static_cast< DiscreteSatelliteBlock * >( v_Block[ i+1 ] )->i2p_y( o , l ), -1.0 ));
      symmetry[i].set_function( new LinearFunction( std::move( v_vars )));
      symmetry[i].set_rhs( 0.0 );
      symmetry[i].set_lhs( -Inf< double >() );
      }
   }

   //////add_static_constraint( symmetry , "symmetry" );

   AR = true;
   std::cout << "Constraints charged!\n";
}

 }  // end( DiscreteConstellationBlock::generate_abstract_constraints() )


/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR CHECKING THE DiscreteConstellationBlock -----*/
/*--------------------------------------------------------------------------*/
// checks feasibility of the DiscreteConstellationBlock-level constraints
// only (thetaM, observation, observation1), cf. ConstellationBlock::is_feasible()

bool DiscreteConstellationBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 1e-1;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
  }
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return( 
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  RowConstraint::is_feasible( thetaM , tol , rel_viol )
  && RowConstraint::is_feasible( observation , tol , rel_viol )
  && RowConstraint::is_feasible( observation1, tol , rel_viol )
);

}  // end( DiscreteConstellationBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
// TODO: printing the DiscreteConstellationBlock instance is not implemented yet

void DiscreteConstellationBlock::print( std::ostream & output , char vlvl ) const
{
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteConstellationBlock::guts_of_destructor( void )
{
 /* clear() all Constraint to ensure that they do not bother to un-register
    themselves from Variable that are going to be deleted anyway. Then
    deletes all the "abstract representation", if any. */

 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 }  // end( guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*---------------------- End File DiscreteConstellationBlock.cpp -----------*/
/*--------------------------------------------------------------------------*/
