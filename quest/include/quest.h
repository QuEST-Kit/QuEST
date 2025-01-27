/** @file
 * The main QuEST header, exposing the entire API.
 */

#ifndef QUEST_H
#define QUEST_H


#include "quest/include/calculations.h"
#include "quest/include/debug.h"
#include "quest/include/decoherence.h"
#include "quest/include/environment.h"
#include "quest/include/initialisations.h"
#include "quest/include/channels.h"
#include "quest/include/modes.h"
#include "quest/include/operations.h"
#include "quest/include/paulis.h"
#include "quest/include/precision.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/types.h"
#include "quest/include/version.h"
#include "quest/include/wrappers.h"


#if INCLUDE_DEPRECATED_FUNCTIONS
    #include "quest/include/deprecated.h"
#endif



#endif // QUEST_H