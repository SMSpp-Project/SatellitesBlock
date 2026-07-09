/*--------------------------------------------------------------------------*/
/*---------------------------- File DiscreteConstellationBlock.h -----------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteConstellationBlock, which implements the
 * Block concept [see Block.h] for a Multicommodity Min Cost Flow problem.
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

#ifndef __DiscreteConstellationBlock
 #define __DiscreteConstellationBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "DiscreteSatelliteBlock.h"
#include "ColVariable.h"
#include "FRowConstraint.h"
#include "Configuration.h"
#include "Objective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*----------------------- DiscreteConstellationBlock-RELATED TYPES ---------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from DiscreteSatelliteBlock.
 *
 *  @{ */

 using CNumber = DiscreteSatelliteBlock::CNumber;
 using c_RHSValue = RowConstraint::c_RHSValue;
 using Vec_CNumber = DiscreteSatelliteBlock::Vec_CNumber;
 using FNumber = DiscreteSatelliteBlock::FNumber;
 using Vec_FNumber = DiscreteSatelliteBlock::Vec_FNumber;

 using FMultiVector = std::vector< Vec_FNumber >;
 using CMultiVector = std::vector< Vec_CNumber >;
 using MultiSubset = std::vector< Block::Subset >;

 using Vec_Bool = std::vector< bool >;

/** @}  end( types ) */
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup DiscreteConstellationBlock_CLASSES Classes in DiscreteConstellationBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS DiscreteConstellationBlock --------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple MMCF Block concept.

class DiscreteConstellationBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 enum MCFType { kMCF , kSPT };

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of DiscreteConstellationBlock
 /** Constructor of DiscreteConstellationBlock. It accepts a pointer to the father
  * Block, which can be of any type. */

 DiscreteConstellationBlock( Block *father = nullptr ) : Block( father ) , AR( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of DiscreteConstellationBlock

 virtual ~DiscreteConstellationBlock() { guts_of_destructor(); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance from the given file in the given format

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the DiscreteConstellationBlock out of an istream

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the "abstract representation" of the Variable of the Block

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/**@} ----------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE DiscreteConstellationBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the DiscreteConstellationBlock
 *  @{ */

 /// print the DiscreteConstellationBlock on an ostream with the given verbosity
 /** Print the DiscreteConstellationBlock on an ostream. So far vlvl is ignored and only very
  * basic information is printed.
  *
  * TODO: implement some verbosity level that produce output files in at
  *       least some of the single-file formats supported by load(); note
  *       that for multi-file formats, print( std::string & ) must be used.
  */

 void print( std::ostream & output , char vlvl = 0 ) const override;
  
/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the DiscreteSatelliteBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the DiscreteConstellationBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( f_sense );
  }

  Index get_numSat( void ) const { return( satellites ); }

  Index get_orbits( int j ) const { return( static_cast< DiscreteSatelliteBlock * >( v_Block[ j ] )->get_orbits() ); }
  
  Index get_ell( int j ) const { return( static_cast< DiscreteSatelliteBlock * >( v_Block[ j ] )->get_ell() ); }
  
  FNumber get_solution( int k, int j, int tt ) const { 
    return((static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->i2p_y(j,tt))->get_value() ); 
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
 *  @{ */
/* @} ----------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 ///< the static strong forcing constrs
 
 int f_sense = Objective::eMin;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );
 
/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert it in the Block factory

/*--------------------------------------------------------------------------*/

 FNumber satellites;   /// the maximum number of satellites in the constellation
 FNumber targets;      ///< the number of targets
 FNumber time_step;    ///< the time discretization step
 FNumber horizon;      ///< the simulation horizon [sec]
  
 Vec_CNumber indexOrbitSat;  
 Vec_CNumber ellSat; 
 Vec_CNumber periods;      ///< the number of periods per target

 Vec_FNumber thetaValF;           
 Vec_FNumber ThetaValMaxSat;
 Vec_CNumber altitude;

 boost::multi_array< double , 5 > obs;
 boost::multi_array< double , 2 > thetaValSat; 

 FRowConstraint thetaM;  
  
 std::vector< FRowConstraint > symmetry; 

 boost::multi_array< FRowConstraint , 2 > observation; /// the observation constraints
 boost::multi_array< FRowConstraint , 2 > observation1; /// the observation1 constraints

 bool AR;
  
 };  // end( class( DiscreteConstellationBlock ) )

/*--------------------------------------------------------------------------*/

/*@}  end( group( DiscreteConstellationBlock_CLASSES ) ) -------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* DiscreteConstellationBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File DiscreteConstellationBlock.h -------------*/
/*--------------------------------------------------------------------------*/

