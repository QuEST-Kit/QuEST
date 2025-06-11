/** @file
 * The main QuEST header, exposing the entire API.
 * This header is intendedly included by user
 * source-code, and is both C11 and C++14 compatible.
 * Preprocessor 'INCLUDE_DEPRECATED_FUNCTIONS' can
 * be defined as 1 to additionally include QuEST's
 * deprecated v3 API, before including this header.
 * 
 * @author Tyson Jones
 * @author Luc Jaulmes (patching CMake install)
 * 
 * @defgroup api 📋 API
 */

/**
 * @page apilink 📋 API
 * The API documentation can be viewed at @ref api.
 * 
 * We're working hard to move that page up one level. 😎
 */

/**
 * @page testlink 🧪 Tests
 * 
 * The unit and integration tests can be viewed at @ref tests.
 * 
 * We're working hard to move that page up one level. 😎
 */

#ifndef QUEST_H
#define QUEST_H


// bind compile settings to installed exec
#if !0
#define FLOAT_PRECISION 2
#define COMPILE_MPI 0
#define COMPILE_OPENMP 0
#define COMPILE_CUDA 0
#define COMPILE_CUQUANTUM 0
#endif


// include version first so it is accessible to 
// debuggers in case a subsequent include fails
#include "quest/include/version.h"

// include before API headers since it validates
// preprocessor configuration, and affirms macro
// preconditions assumed by subsequent header
#include "quest/include/modes.h"

#include "quest/include/precision.h"
#include "quest/include/types.h"
#include "quest/include/calculations.h"
#include "quest/include/debug.h"
#include "quest/include/decoherence.h"
#include "quest/include/environment.h"
#include "quest/include/initialisations.h"
#include "quest/include/channels.h"
#include "quest/include/operations.h"
#include "quest/include/paulis.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/wrappers.h"


#if INCLUDE_DEPRECATED_FUNCTIONS
    #include "quest/include/deprecated.h"
#endif



#endif // QUEST_H
