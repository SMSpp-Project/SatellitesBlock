/*--------------------------------------------------------------------------*/
/*---------------------------- File MultiTargetBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiTargetBlock, which implements the
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

#ifndef __MultiTargetBlock
 #define __MultiTargetBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "SingleTargetBlock.h"
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
/*----------------------- MultiTargetBlock-RELATED TYPES --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from SingleTargetBlock.
 *
 *  @{ */

 using CNumber = SingleTargetBlock::CNumber;
 using c_RHSValue = RowConstraint::c_RHSValue;
 using Vec_CNumber = SingleTargetBlock::Vec_CNumber;
 using FNumber = SingleTargetBlock::FNumber;
 using Vec_FNumber = SingleTargetBlock::Vec_FNumber;

 using FMultiVector = std::vector< Vec_FNumber >;
 using CMultiVector = std::vector< Vec_CNumber >;
 using MultiSubset = std::vector< Block::Subset >;

 using Vec_Bool = std::vector< bool >;

/** @}  end( types ) */
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MultiTargetBlock_CLASSES Classes in MultiTargetBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MultiTargetBlock ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple MMCF Block concept.

class MultiTargetBlock : public Block
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

 /// constructor of MultiTargetBlock
 /** Constructor of MultiTargetBlock. It accepts a pointer to the father
  * Block, which can be of any type. */

 MultiTargetBlock( Block *father = nullptr ) : Block( father ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of MultiTargetBlock

 virtual ~MultiTargetBlock() { guts_of_destructor(); }

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
  * If there is any Solver attached to this MultiTargetBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the MultiTargetBlock out of an istream
 /** Load the MultiTargetBlock out of an istream. Handles the two single-file
  * formats, i.e., Canad and PPRN.
  *
  * TODO: properly document the formats.
  *
  * If there is any Solver attached to this MultiTargetBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the "abstract representation" of the Variable of the Block
 /** This method generates the "abstract representation" of the Variable of
  * the MultiTargetBlock, and in fact it decides which formulation of the MMCF
  * problem is implemented. This is controlled by the parameter stvv. If stvv
  * is not nullptr and it is a SimpleConfiguration< int >, or if
  * f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
  * SimpleConfiguration< int >, then the f_value (an int) dictates which
  * MMCF formulation as follows:
  *
  * - [1]: the standard knapsack formulation in which get_NArcs()
  *   BinaryKnapsackBlock sub-Block are constructed, one for each commodity,
  *   and the flow constraints are handled in the father MultiTargetBlock;
  *
  * - [0]: the standard flow formulation in which get_NComm() SingleTargetBlock
  *   sub-Block are constructed, one for each commodity, and the
  *   linking constraints are handled in the father MultiTargetBlock;
  *
  * - [other ones possibly to follow].
  * 
  *  by default is considered the Flow relaxation
  */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE MultiTargetBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the MultiTargetBlock
 *  @{ */

 /// print the MultiTargetBlock on an ostream with the given verbosity
 /** Print the MultiTargetBlock on an ostream. So far vlvl is ignored and only very
  * basic information is printed.
  *
  * TODO: implement some verbosity level that produce output files in at
  *       least some of the single-file formats supported by load(); note
  *       that for multi-file formats, print( std::string & ) must be used.
  */

 void print( std::ostream & output , char vlvl = 0 ) const override;
  
/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SingleTargetBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the MultiTargetBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// getting the current sense of the Objective, which is minimization

 int get_objective_sense( void ) const override final {
  return( f_sense );
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

 FNumber satellites;   ///< the number of satellites in the MultiTarget
 FNumber targets;      ///< the number of targets
 FNumber time_step;    ///< the time discretization step
 FNumber horizon;      ///< the simulation horizon [sec]
 FNumber indexOrbit;      ///< the number of possible orbital plane

 boost::multi_array< FRowConstraint , 3 > duplicate_pi; ///< duplicate_pi constraints; 
  
 };  // end( class( MultiTargetBlock ) )

/*--------------------------------------------------------------------------*/

/*@}  end( group( MultiTargetBlock_CLASSES ) ) ---------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MultiTargetBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiTargetBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/

