/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to call to effect local simulation subroutines 
 * upon Quregs. The data needed for these subroutines must 
 * already be localised into the appropriate memory (RAM vs VRAM)
 * and location (qureg's amplitudes or buffer space), as is
 * performed by localiser.cpp. These subroutines are ergo
 * embarrassingly parallel.
 */

#ifndef ACCELERATOR_HPP
#define ACCELERATOR_HPP



#endif // ACCELERATOR_HPP