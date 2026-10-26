/*--------------------------------------------------------------------------*/
/*----------------------- File SatellitesBlock.cppm ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.satellitesblock, which exports every
 * name that the headers of SatellitesBlock declare at namespace scope,
 * together with the modules it depends on; a translation unit can then write
 *
 *     import smspp.satellitesblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "ConstellationBlock.h"
#include "DiscreteConstellationBlock.h"
#include "DiscreteSatelliteSolver.h"
#include "MultiTargetBlock.h"
#include "MultiTargetBlockv2.h"
#include "SatelliteSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.satellitesblock;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // ConstellationBlock.h
 using SMSpp_di_unipi_it::CMultiVector;
 using SMSpp_di_unipi_it::CNumber;
 using SMSpp_di_unipi_it::ConstellationBlock;
 using SMSpp_di_unipi_it::FMultiVector;
 using SMSpp_di_unipi_it::FNumber;
 using SMSpp_di_unipi_it::MultiSubset;
 using SMSpp_di_unipi_it::Vec_Bool;
 using SMSpp_di_unipi_it::Vec_CNumber;
 using SMSpp_di_unipi_it::Vec_FNumber;
 using SMSpp_di_unipi_it::c_RHSValue;

 // DiscreteConstellationBlock.h
 using SMSpp_di_unipi_it::DiscreteConstellationBlock;

 // DiscreteSatelliteBlock.h
 using SMSpp_di_unipi_it::DiscreteSatelliteBlock;
 using SMSpp_di_unipi_it::DiscreteSatelliteSolution;
 using SMSpp_di_unipi_it::Vec_DiscreteSatelliteBlock;
 using SMSpp_di_unipi_it::Vec_DiscreteSatelliteBlock_it;
 using SMSpp_di_unipi_it::c_Vec_DiscreteSatelliteBlock;
 using SMSpp_di_unipi_it::c_Vec_DiscreteSatelliteBlock_it;
 using SMSpp_di_unipi_it::p_DiscreteSatelliteBlock;

 // DiscreteSatelliteSolver.h
 using SMSpp_di_unipi_it::DiscreteSatelliteSolver;

 // MultiTargetBlock.h
 using SMSpp_di_unipi_it::MultiTargetBlock;

 // MultiTargetBlockv2.h
 using SMSpp_di_unipi_it::MultiTargetBlockv2;
 using SMSpp_di_unipi_it::MultiTargetSolution;
 using SMSpp_di_unipi_it::Vec_MultiTargetBlockv2;
 using SMSpp_di_unipi_it::Vec_MultiTargetBlockv2_it;
 using SMSpp_di_unipi_it::c_Vec_MultiTargetBlockv2;
 using SMSpp_di_unipi_it::c_Vec_MultiTargetBlockv2_it;
 using SMSpp_di_unipi_it::p_MultiTargetBlockv2;

 // SatelliteBlock.h
 using SMSpp_di_unipi_it::SatelliteBlock;
 using SMSpp_di_unipi_it::SatelliteSolution;
 using SMSpp_di_unipi_it::Vec_SatelliteBlock;
 using SMSpp_di_unipi_it::Vec_SatelliteBlock_it;
 using SMSpp_di_unipi_it::c_Vec_SatelliteBlock;
 using SMSpp_di_unipi_it::c_Vec_SatelliteBlock_it;
 using SMSpp_di_unipi_it::p_SatelliteBlock;

 // SatelliteSolver.h
 using SMSpp_di_unipi_it::SatelliteSolver;

 // SingleTargetBlock.h
 using SMSpp_di_unipi_it::SingleTargetBlock;
 using SMSpp_di_unipi_it::SingleTargetSolution;
 using SMSpp_di_unipi_it::Vec_SingleTargetBlock;
 using SMSpp_di_unipi_it::Vec_SingleTargetBlock_it;
 using SMSpp_di_unipi_it::c_Vec_SingleTargetBlock;
 using SMSpp_di_unipi_it::c_Vec_SingleTargetBlock_it;
 using SMSpp_di_unipi_it::p_SingleTargetBlock;
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File SatellitesBlock.cppm ----------------------*/
/*--------------------------------------------------------------------------*/
