/*--------------------------------------------------------------------------*/
/*---------------------------- File ConstellationBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class ConstellationBlock, which implements the
 * Block concept [see Block.h] for Satellite Constellation Design Problem (SCDP).
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
/*-------------------------- CLASS ConstellationBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of a simple ConstellationBlock concept.
/* ConstellationBlock is composed of some SatelliteBlocks linked by the observability 
constraints for the targets, indicating that each target has to be observed by at
least one satellite of the constellation within the revisit time 
* \f[
    \sum_{\substack{t \in T(k,\Delta t_m,dt) \\ i \in \mathcal [s]}} \xi_{i,t,m} \ge 1 \quad 
    \forall k \in [\lfloor T\slash\Delta t_m\rfloor], \forall m \in \mathcal{X},
 * \f]
where \Delta t_m is the revisit time for target m \in \mathcal{X} and [s] = \{1,2,...,s\} is 
the set of the satellites in the constellation. Variables \xi_{i,t,m} indicating whether
the sallite i observes target m at time stamp t (see SatelliteBlock).

Moreover, we implement a family of symmetry breaking constraints, i.e, z_{i} \leq z_{i+1},
meaning that, if the (i+1)-th satellite is active in the constellation, then also 
the previous satellites should be active in the constellation. These latter constraints
have been implemented here since we consider the case where all the satellites are 
identical, i.e., they are caracterised by the same target capturing technology.
*/

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

 ConstellationBlock( Block *father = nullptr ) : Block( father ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of ConstellationBlock

 virtual ~ConstellationBlock() { guts_of_destructor(); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the instance from the given file in the given format
 /** Loads a SDCP instance using filename as the "base filename".
  * If there is any Solver attached to this ConstellationBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( const std::string & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// load the ConstellationBlock out of an istream
 /*
  * If there is any Solver attached to this ConstellationBlock then a NBModification
  * (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// generate the "abstract representation" of the Variable of the Block
 /** This method generates the "abstract representation" of the Variable of
  * the ConstellationBlock. 
  */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE ConstellationBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the ConstellationBlock
 *  @{ */

 /// print the ConstellationBlock on an ostream with the given verbosity

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
 FNumber time_step;    ///< the time discretization step [seconds]
 FNumber horizon;      ///< the simulation horizon [seconds]

 Vec_CNumber periods;      ///< the number of periods per target

 boost::multi_array< FRowConstraint , 2 > observation; /// the observation constraints
 FRowConstraint observation1; /// the observation1 constraints
 FRowConstraint observation2; /// the observation2 constraints

 std::vector< FRowConstraint > symmetry; 
  
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

