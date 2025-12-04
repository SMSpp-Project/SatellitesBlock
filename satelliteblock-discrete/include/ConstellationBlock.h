/*--------------------------------------------------------------------------*/
/*---------------------------- File ConstellationBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class ConstellationBlock, which implements the
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

#ifndef __ConstellationBlock
 #define __ConstellationBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "SatelliteBlock.h"
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
/*----------------------- ConstellationBlock-RELATED TYPES --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from SatelliteBlock.
 *
 *  @{ */

 using CNumber = SatelliteBlock::CNumber;
 using c_RHSValue = RowConstraint::c_RHSValue;
 using Vec_CNumber = SatelliteBlock::Vec_CNumber;
 using FNumber = SatelliteBlock::FNumber;
 using Vec_FNumber = SatelliteBlock::Vec_FNumber;

 using FMultiVector = std::vector< Vec_FNumber >;
 using CMultiVector = std::vector< Vec_CNumber >;
 using MultiSubset = std::vector< Block::Subset >;

 using Vec_Bool = std::vector< bool >;

/** @}  end( types ) */
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ConstellationBlock_CLASSES Classes in ConstellationBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ConstellationBlock ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple MMCF Block concept.

class ConstellationBlock : public Block
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

 /// constructor of ConstellationBlock
 /** Constructor of ConstellationBlock. It accepts a pointer to the father
  * Block, which can be of any type. */

 ConstellationBlock( Block *father = nullptr ) : Block( father ) , AR( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of ConstellationBlock

 virtual ~ConstellationBlock() { guts_of_destructor(); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance from the given file in the given format
 /** Loads a MMCF instance using filename as the "base filename". This method
  * supports several formats depending on \p frmt, that is case-insensitive.
  * In particular, for two single-file formats
  *
  * - frmt == 0 (default) or frmt == 'c': PPRN format
  *
  * - frmt == 's': Canad format
  *
  * it behaves just as the Block method (just open an ifstream and
  * dispatch it load( std::istream & ). However, it also supports 5
  * multi-file formats:
  *
  * - 'm': Mnetgen format
  * - 'p': Jones-Lustig PSP (product-specific problem) format
  * - 'o': Jones-Lustig OSP (origin-specific problem) format
  * - 'd': Jones-Lustig OSP (origin-destination problem) format
  * - 'u': same as 'd' but supply information is looked at in file
  *        input + ".od" rather than input + ".sup" as in all the
  *        other cases
  *
  * where input (prefixed as set by set_filename_prefix(), if any) is
  * completed by the appropriate suffixes ".nod", ".arc", ".mut", ".sup"
  * or ".od" to load different parts of the description of the MMCF
  * instance.
  *
  * TODO: properly document all the formats.
  *
  * If there is any Solver attached to this ConstellationBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the ConstellationBlock out of an istream
 /** Load the ConstellationBlock out of an istream. Handles the two single-file
  * formats, i.e., Canad and PPRN.
  *
  * TODO: properly document the formats.
  *
  * If there is any Solver attached to this ConstellationBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the "abstract representation" of the Variable of the Block
 /** This method generates the "abstract representation" of the Variable of
  * the ConstellationBlock, and in fact it decides which formulation of the MMCF
  * problem is implemented. This is controlled by the parameter stvv. If stvv
  * is not nullptr and it is a SimpleConfiguration< int >, or if
  * f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
  * SimpleConfiguration< int >, then the f_value (an int) dictates which
  * MMCF formulation as follows:
  *
  * - [1]: the standard knapsack formulation in which get_NArcs()
  *   BinaryKnapsackBlock sub-Block are constructed, one for each commodity,
  *   and the flow constraints are handled in the father ConstellationBlock;
  *
  * - [0]: the standard flow formulation in which get_NComm() SatelliteBlock
  *   sub-Block are constructed, one for each commodity, and the
  *   linking constraints are handled in the father ConstellationBlock;
  *
  * - [other ones possibly to follow].
  * 
  *  by default is considered the Flow relaxation
  */

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
/*--------------- METHODS FOR PRINTING & SAVING THE ConstellationBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the ConstellationBlock
 *  @{ */

 /// print the ConstellationBlock on an ostream with the given verbosity
 /** Print the ConstellationBlock on an ostream. So far vlvl is ignored and only very
  * basic information is printed.
  *
  * TODO: implement some verbosity level that produce output files in at
  *       least some of the single-file formats supported by load(); note
  *       that for multi-file formats, print( std::string & ) must be used.
  */

 void print( std::ostream & output , char vlvl = 0 ) const override;
  
/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SatelliteBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the ConstellationBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( f_sense );
  }

  Index get_numSat( void ) const { return( satellites ); }

  Index get_orbits( int j ) const { return( static_cast< SatelliteBlock * >( v_Block[ j ] )->get_orbits() ); }
  Index get_ell( int j ) const { return( static_cast< SatelliteBlock * >( v_Block[ j ] )->get_ell() ); }
  FNumber get_solution( int k, int j, int tt ) const { return((static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_y(j,tt))->get_value() ); }

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
 *  @{ */


/*--------------------------------------------------------------------------*/
 /** called at the end of any constructor, does some initializations that are
  * common to them all: it is "protected" for allowing derived classes that
  * use the "void" constructor to call it. */

 //void CmnIntlz( void );

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
  
 };  // end( class( ConstellationBlock ) )

/*--------------------------------------------------------------------------*/

/*@}  end( group( ConstellationBlock_CLASSES ) ) ---------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ConstellationBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ConstellationBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/

