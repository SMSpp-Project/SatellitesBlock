/*--------------------------------------------------------------------------*/
/*------------------------- File MultiTargetBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MultiTargetBlock class.
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

#include "MultiTargetBlock.h"

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

// register MultiTargetBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( MultiTargetBlock );

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

void  MultiTargetBlock::load( std::istream & input , char frmt ){
}

void MultiTargetBlock::load( const std::string & input , char frmt )
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
 Vec_CNumber altitude;
 Vec_CNumber periodSat;
 int indexLen = 0;
 double j = 0.0;

 for( Index i = 0 ; i < horizon/3600.0 ; ++i ){
   j++;
	altitudeSetVal = cbrt( ( MU* pow((horizon)/j,2.0) ) / (4.0*pow(PI,2.0))) - RAYON;
	if ( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 )
   {
		indexLen++;
   }
 }

 FNumber altSet = indexLen;
 altitude.resize( indexLen );
 periodSat.resize( indexLen );
 indexLen = 0;

 j = 0.0;
 for( Index i = 0 ; i < horizon/3600.0 ; ++i )
 {
   j++;
	altitudeSetVal = cbrt( ( MU * pow((horizon)/j, 2.0) ) / (4.0*pow(PI, 2.0) )) - RAYON;
   //std::cout << altitudeSetVal << "\n";
	if ( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 )
   {
		altitude[indexLen] = altitudeSetVal;
		periodSat[indexLen] = horizon/j;
      indexLen++;
   }
 }

 Latitude.resize( targets );
 Longitude.resize( targets );
 for( Index i = 0 ; i < targets ; ++i )
 {
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

 int alpha_lim_flag = 0;
 
 for( Index ii = 0 ; ii < altSet ; ++ii ){
   if(((RAYON+altitude[ ii ])/RAYON)*sin(aHalf) > 1){
      aHalf = asin((RAYON/(RAYON+altitude[ ii ])));
      std::cout << "WARNING: computed alpha_lim\n";
      alpha_lim_flag = 1;
      break;
   } 
 }

 double alt = altitude[altSet-1];
 //double alt = altitude[0];
 for( Index ii = 0 ; ii < altSet ; ++ii ){
   t_p[ ii ] = sqrt(MU / ( pow(RAYON+altitude[ii],3)));
   t_u[ ii ] = sqrt(MU / (RAYON+altitude[ ii ]));
   t_GM[ ii ] = sqrt((RAYON+altitude[ ii ]) / MU);
   thetaVal[ ii ] = -aHalf + asin(((RAYON+altitude[ii])/RAYON)*sin(aHalf));
 }

 //double Theta_max = thetaVal[altSet-1];
 std::cout << "alpha_flag=" << alpha_lim_flag << "\n";
 
 /***********************/
 if (alpha_lim_flag == 1) {
   /////Theta_min = thetaVal[altSet-1];//*1.2;
   std::cout << "theta_min: " << Theta_min << "\n";
 }
 /***********************/
 
 //Theta_max = Theta_min;
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
 //double altitudeFinal = altitude[0];
 
 //std::cout << "theta_min: " << Theta_min << "\n";
 std::cout << "alpha_half: " << aHalf*180/PI << "\n";
 std::cout << "theta: " << thetaValFinal << "\n";
 std::cout << "altitude: " << altitudeFinal << "\n";

 FNumber t = horizon / time_step;

 boost::multi_array< double , 3 > CoverageSatLat(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLat1(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong1(boost::extents[targets][t][incSet*ascSet*anmSet]);

 double lat_Sat;
 double long_Sat;

 //for( Index ii = 0 ; ii < altSet ; ++ii )
 //{
 int ii = indexLen-1;
 //int ii = 0; 
 int index1 = -1;
 indexOrbit = 0;
 int addOrbit = 0;
 int indexOrbit1 = 0;
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
                CoverageSatLat1[ i ][ j ][ index1 ] = std::abs(Latitude[i] - lat_Sat); 
                CoverageSatLong1[ i ][ j ][ index1 ] = std::abs(Longitude[i] - long_Sat) * cos(Latitude[i]); 
                if(cos(Latitude[i]) < 0)
                  std::cout << "ERROR!" << "\n";
                //std::cout << i << " " << j << " " << ii << " " << jj << " " << k << " " << l << " " << "\n";
                if (CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValFinal and CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValFinal)
                  indexOrbit1 += 1;
             }
            }
            if (indexOrbit1>=1){
               for( Index j = 0 ; j < t ; ++j ) 
               {
                  for( Index i = 0 ; i < targets ; ++i ) 
                  {
                  CoverageSatLat[ i ][ j ][ indexOrbit ] = CoverageSatLat1[ i ][ j ][ index1 ];
                  CoverageSatLong[ i ][ j ][ indexOrbit ] = CoverageSatLong1[ i ][ j ][ index1 ];
                  }
               }
               indexOrbit += 1;
            }
            
          }
       }
    }
 //}

 iFile >> satellites;
 std::cout << "number Orbits: " << indexOrbit << "\n";

 v_Block.resize( targets );
 std::cout << "SATELLITES: " << satellites << "\n";
  
 for( Index i = 0 ; i < targets ; ++i ){
   boost::multi_array< double , 2 > CoverageSatLatThis(boost::extents[t][indexOrbit]);
   boost::multi_array< double , 2 > CoverageSatLongThis(boost::extents[t][indexOrbit]);
   
   for( Index j = 0 ; j < t ; ++j )  {
      for( Index k = 0 ; k < indexOrbit ; ++k ) { 
         CoverageSatLatThis[ j ][ k ] = CoverageSatLat[ i ][ j ][ k ];
         CoverageSatLongThis[ j ][ k ] = CoverageSatLong[ i ][ j ][ k ];
      }
   }

    auto SB = new SingleTargetBlock( this );
    SB->load( satellites , time_step , horizon , altitudeFinal , thetaValFinal ,
               indexOrbit , aHalf , CoverageSatLatThis , CoverageSatLongThis );
    v_Block[ i ] = SB;
 }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( MultiTargetBlock::load( const std::string & input )

/*--------------------------------------------------------------------------*/

void MultiTargetBlock::generate_abstract_variables( Configuration * stvv )
{

  for( auto blck : v_Block )
    blck->generate_abstract_variables();

  std::cout << "Variables generated!\n";

 }

/*--------------------------------------------------------------------------*/

void MultiTargetBlock::generate_abstract_constraints( Configuration * stcc )
{

   for( auto blck : v_Block )
    blck->generate_abstract_constraints();

   duplicate_pi.resize( boost::multi_array_types::extent_gen()[ targets - 1  ][ satellites ][ indexOrbit ] );
   for( Index i = 0 ; i < targets - 1  ; ++i ) {
      for( Index j = 0 ; j < satellites  ; ++j ) {
         for( Index k = 0 ; k < indexOrbit  ; ++k ) {
            LinearFunction::v_coeff_pair v_vars;
            v_vars.push_back( std::make_pair( static_cast< SingleTargetBlock * >( v_Block[ i ] )->i2p_pi( j , k ), 1.0 ));
            v_vars.push_back( std::make_pair( static_cast< SingleTargetBlock * >( v_Block[ i+1 ] )->i2p_pi( j , k ), -1.0 ));
            duplicate_pi[i][j][k].set_function( new LinearFunction( std::move( v_vars )));
            duplicate_pi[i][j][k].set_rhs( 0.0 ); 
            duplicate_pi[i][j][k].set_lhs( 0.0 );
         }
      }
   }

   add_static_constraint( duplicate_pi , "duplicate_pi" );

   std::cout << "Constraints generated!\n";

  // generate the observability constraints  - - - - - - - - - - - - - - -

 }  // end( MultiTargetBlock::generate_abstract_constraints() )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlock::print( std::ostream & output , char vlvl ) const
{
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlock::guts_of_destructor( void )
{
 /* clear() all Constraint to ensure that they do not bother to un-register
    themselves from Variable that are going to be deleted anyway. Then
    deletes all the "abstract representation", if any. */

 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 // not needed, there isn't any - reset_static_variables();
 // not needed, there isn't any - reset_dynamic_constraints();
 // not needed, there isn't any - reset_dynamic_variables();
 // not needed, there isn't any - reset_objective();

 }  // end( guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiTargetBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
