/** @file
 * The main QuEST header, exposing the entire API.
 * This header is intendedly included by user
 * source-code, and is both C11 and C++14 compatible.
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

// include config.h first to define macros
// consulted by subsequent headers
#include "quest/include/config.h"

#include "quest/include/modes.h"
#include "quest/include/precision.h"
#include "quest/include/types.h"
#include "quest/include/calculations.h"
#include "quest/include/debug.h"
#include "quest/include/decoherence.h"
#include "quest/include/environment.h"
#include "quest/include/trotterisation.h"
#include "quest/include/initialisations.h"
#include "quest/include/channels.h"
#include "quest/include/multiplication.h"
#include "quest/include/operations.h"
#include "quest/include/paulis.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/wrappers.h"


#if INCLUDE_DEPRECATED_FUNCTIONS
    #include "quest/include/deprecated.h"
#endif



#endif // QUEST_H
