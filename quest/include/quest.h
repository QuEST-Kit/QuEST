/** @file
 * The main QuEST header, exposing the entire API.
 * This header is intendedly included by user
 * source-code, and is both C11 and C++14 compatible.
 * Preprocessor 'INCLUDE_DEPRECATED_FUNCTIONS' can
 * be defined as 1 to additionally include QuEST's
 * deprecated v3 API, before including this header.
 * 
 * @author Tyson Jones
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


// include version first so it is accessible to 
// debuggers in case a subsequent include fails
#include "quest/version.h"

// include before API headers since it validates
// preprocessor configuration, and affirms macro
// preconditions assumed by subsequent header
#include "quest/modes.h"

#include "quest/precision.h"
#include "quest/types.h"
#include "quest/calculations.h"
#include "quest/debug.h"
#include "quest/decoherence.h"
#include "quest/environment.h"
#include "quest/initialisations.h"
#include "quest/channels.h"
#include "quest/operations.h"
#include "quest/paulis.h"
#include "quest/qureg.h"
#include "quest/matrices.h"
#include "quest/wrappers.h"


#if INCLUDE_DEPRECATED_FUNCTIONS
    #include "quest/deprecated.h"
#endif



#endif // QUEST_H