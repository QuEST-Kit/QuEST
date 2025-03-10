/** @file
 * Functions which unambiguously identify the buffer indices
 * at which communicated amplitudes are sent and received.
 * These are inlined mainly to avoid symbol duplication as
 * a header-only file, but also so that callers of e.g.
 * getBufferRecvInd() (which use the result as an index-offset 
 * in hot loops) can exploit the compile-time known constant.
 * 
 * @author Tyson Jones
 */

#ifndef COMM_INDICES_HPP
#define COMM_INDICES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include <utility>



static inline qindex getSubBufferSendInd(Qureg qureg) {

    // the maximum size of a swapped sub-buffer is half its capacity, so we
    // will always pack sub-buffers starting from half capacity
    return qureg.numAmpsPerNode / 2;
}


constexpr static inline qindex getBufferRecvInd() {

    // we always receive amplitudes to the start of the buffer, regardless
    // of whether we are receiving a full or sub-buffer
    return 0;
}


static inline std::pair<qindex,qindex> getSubBufferSendRecvInds(Qureg qureg) {

    return {getSubBufferSendInd(qureg), getBufferRecvInd()};
}



#endif // COMM_INDICES_HPP