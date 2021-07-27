/*
 * =====================================================================================
 *
 *       Filename:  Interface.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Saturday 13 July 2019 12:51:38  PDT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#if brgr==1
#include "../Burger/Mean.hpp"
#endif
#if quasi_geo==1
#include "../QG/Mean.hpp"
#endif
#if need_locaInterface == 1
#include "NOX_Epetra_LinearSystem_Amesos.H"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "../swe/FVM_LocaInterface.H"
#include "../swe/ImplicitTimeStepper.H"
#include "../swe/ThetaStepperEvaluator.H"
#endif
