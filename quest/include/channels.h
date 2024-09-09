/*
 * STRUCTS
 */


typedef struct {

    // fields are const to prevent user modification
    const int numQubits;
    const qindex numRows;
    
    qcomp** cpuElems;
    qcomp* gpuElemsFlat;

    // whether the user has ever synchronised memory to the GPU, which is performed automatically
    // when calling functions like setCompMatr(), but which requires manual invocation with
    // syncCompMatr() after manual modification of the cpuElem. Note this can only indicate whether
    // the matrix has EVER been synced; it cannot be used to detect whether manual modifications
    // made after an initial sync have been re-synched. This is a heap pointer to remain mutable.
    int* const wasGpuSynced;

} SuperOp;


typedef struct {

    // fields are const to prevent user modification
    const int numQubits;

    // representation of the map as a collection of Kraus operators, kept exclusively 
    // in CPU memory, and used only for CPTP validation and reporting the map
    const qindex numMatrices;
    const qindex numRows;
    qcomp*** matrices;

    // representation of the map as a single superoperator, used for simulation
    SuperOp superop;

    // CPTP-ness is determined at validation; 0 or 1, or -1 to indicate unknown. The flag is 
    // stored in heap so even copies of structs are mutable, but pointer itself is immutable.
    int* const isCPTP;

} KrausMap;


// we define no fixed-size versions (e.g. KrausMap1/2), unlike we did for CompMatr1/2
// and DiagMatr1/2. This is because the 2-qubit superoperator is 256 elements big, and
// seems inadvisably large to be passing-by-copy through the QuEST backend layers, and
// would need explicit GPU memory allocation at each invocation of mixKrausMap2() (it
// exceeds the max number of CUDA kernel args). Furthermore, repeatedly calling
// createKrausMap2() would repeatedly invoke ~256*16 flops to compute te superoperator,
// which may be an user-astonishing overhead (more astonishing than the API asymmetry).
// Finally, computing the fixed-size superoperators must be in the header (to avoid
// the issues of qcmop interoperability, just like for getCompMatr1) and could not call
// an inner function which wouldn't be user-exposed; so we would end up redefining the
// superoperator calculation THREE times!



