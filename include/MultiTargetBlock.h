/*--------------------------------------------------------------------------*/
/*---------------------------- File MultiTargetBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiTargetBlock, which implements the
 * Block concept [see Block.h] for a Satellite Constellation Design Problem (SCDP).
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
/*-------------------------- CLASS MultiTargetBlock ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple MultiTargetBlock concept.
/* MultiTargetBlock is composed of some SingleBlocks linked by the 
* configuration constraints for the satellites, indicating that each satellite 
must have a unique orbital configuration
* \f[
    \sum_{i \in [C]} pi[ i ][ c ] = 1, \forall i \in [s],
* \f]
* where [s] = \{1,2,...,s\} is the set of the satellites in the constellation. 
* Variables pi[ i ][ c ] indicating whether the sallite i is in orbital 
configuration c (the possible configuration are pre-computed when loading the 
problem instance, see SingleTargetBlock.h).
*

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
 /** Loads a SCDP instance using filename as the "base filename". 
  * If there is any Solver attached to this MultiTargetBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the MultiTargetBlock out of an istream
 /*
  * If there is any Solver attached to this MultiTargetBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the "abstract representation" of the Variable of the Block
 /** This method generates the "abstract representation" of the Variable of
  * the MultiTargetBlock.
  */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE MultiTargetBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the MultiTargetBlock
 *  @{ */

 /// print the MultiTargetBlock on an ostream with the given verbosity

 void print( std::ostream & output , char vlvl = 0 ) const override;
  
/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the SingleTargetBlock -----*/
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
 FNumber time_step;    ///< the time discretization step [seconds]
 FNumber horizon;      ///< the simulation horizon [seconds]
 FNumber indexOrbit;      ///< the number of possible orbital plane

 //boost::multi_array< FRowConstraint , 3 > duplicate_pi; ///< duplicate_pi constraints; 
 boost::multi_array< FRowConstraint , 2 > duplicate_pi; ///< duplicate_pi constraints; 
  
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

