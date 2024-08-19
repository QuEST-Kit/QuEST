/** @file
 * Functions which unambiguously identify the buffer indices
 * at which communicated amplitudes are sent and received
 */

#ifndef COMM_INDICES_HPP
#define COMM_INDICES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include "quest/src/core/inliner.hpp"

#include <utility>



/*
 * BUFFER INDICES
 *
 * which are inlined mainly to avoid symbol duplication,
 * but also so that callers of getBufferRecvInd() which
 * use the result as an index-offset in hot loops can
 * exploit the compile-time known constant.
 */


INLINE qindex getSubBufferSendInd(Qureg qureg) {

    // the maximum size of a swapped sub-buffer is half its capacity, so we
    // will always pack sub-buffers starting from half capacity
    return qureg.numAmpsPerNode / 2;
}


INLINE qindex getBufferRecvInd() {

    // we always receive amplitudes to the start of the buffer, regardless
    // of whether we are receiving a full or sub-buffer
    return 0;
}


INLINE std::pair<qindex,qindex> getSubBufferSendRecvInds(Qureg qureg) {

    return {getSubBufferSendInd(qureg), getBufferRecvInd()};
}



#endif // COMM_INDICES_HPP