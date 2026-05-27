/** @file
 * A stack-based list of length <= 64, primarily
 * for storing qubit indices, as an alternative to
 * std::vector and associated heap-alloc/copy
 * overheads. Use of List64 optimises few-qubit
 * simulation where STL container costs dominate;
 * and in the GPU backend, use of List64 avoids
 * CUDA memory writes before kernel launches!
 * 
 * This header also defines ConstList64, which is
 * merely 'const List64&', to avoid superfluous
 * stack copies when passing non-mutated List64.
 * 
 * The functions herein are inlined (in this header-
 * only file) in the hopes of unbridled compiler
 * optimisations, but this may prove incompatible
 * with GPU mode (since INLINE specifies __device__,
 * which may be incompatible with initialiser lists)
 * 
 * @author Tyson Jones
 */

#ifndef LISTS_HPP
#define LISTS_HPP

#include "quest/src/core/errors.hpp"
#include "quest/src/core/inliner.hpp"



/*
 * CAPACITY
 *
 * Since stored in stack, we must upperbound the length of
 * a List64; we choose 64, which is around the maximum
 * addressable number of qubits by qindex. In theory, we
 * could permit users to compile-time reduce this length,
 * restricting their max simulable system but speeding up
 * List64 copies in function calls - this may have a
 * measurable benefit for Quregs of 1-8 qubits. But Donald
 * Knuth knows and sees all, and he won't be happy!
 */


constexpr size_t MAX_LIST_LENGTH = 64;



/*
 * LIST64 DECLARATION
 *
 * which mimics an STL container so that it is easily
 * substituted for std::vector in our codebase, but
 * crucially, remains (almost) POD and with no heap
 * allocs, and compatible with CUDA kernels 
 */


struct List64 {

private:

    // Keep data private to dissuade inconsistent
    // access patterns (e.g. .elems vs .data()),
    // and so users cannot invalidly mutate length.
    // Readers may wonder why we avoid std::array;
    // it has a surprise overhead in pass-by-ref!
    int elems[MAX_LIST_LENGTH];

    // We use size_t, over the arguably internally
    // natural int, for consistency with STL containers
    size_t length;

public:

    // Note there is deliberately no constructor!
    // This keeps the struct trivial and compatible
    // with CUDA; we must forego initializer ctors
    // and other syntactic goodies :(

    // let List64 be iterable, e.g. for(auto x : list)
    INLINE auto begin()       { return elems; }
    INLINE auto begin() const { return elems; }
    INLINE auto end()         { return elems + length; }
    INLINE auto end()   const { return elems + length; }

    // let List64 be indexable, e.g. list[3]
    INLINE const int& operator[](int index) const {

        if (index < 0)
            error_smallListIndexWasNegative();
        if (index >= static_cast<int>(length))
            error_smallListIndexExceededLength();

        return elems[index];
    }
    INLINE int& operator[](int index) {

        return const_cast<int&>(
            static_cast<const List64&>(*this)[index]);
    }

    // give List64 all the familiar methods of std::vector
    INLINE void clear() {
        length = 0;
    }
    INLINE bool empty() const { 
        return length == 0; 
    }
    INLINE size_t size() const { 
        return length;
    }
    INLINE int* data() {
        return elems;
    }
    INLINE const int* data() const {
        return elems;
    }

    INLINE void push_back(int elem) {

        if (length >= MAX_LIST_LENGTH)
            error_smallListLengthExceededMax();

        elems[length++] = elem;
    }

    INLINE void resize(size_t newLength, int value=0) {

        if (newLength > MAX_LIST_LENGTH)
            error_smallListLengthExceededMax();

        for (auto i=length; i<newLength; i++)
            elems[i] = value;

        length = newLength;
    }

    INLINE const int& back() const {

        if (length == 0)
            error_smallListWasEmpty();

        return elems[length - 1];
    }
    INLINE int& back() {

        return const_cast<int&>(
            static_cast<const List64&>(*this).back());
    }

    INLINE void assign(size_t count, int value) {

        if (count > MAX_LIST_LENGTH)
            error_smallListLengthExceededMax();

        for (size_t i = 0; i < count; i++)
            elems[i] = value;

        length = count;
    }
};



/*
 * LIST64 CONSTRUCTORS
 *
 * which are separated here because making them actual
 * constructors stops List64 being POD/trivial, and
 * makes it incompatible with CUDA kernels
 */


INLINE List64 lists_getEmptyList64() {

    List64 out{};
    out.clear();
    return out;
}


INLINE List64 lists_getList64(const int* begin, const int* end) {

    if (end < begin)
        error_smallListIndexExceededLength();

    auto length = static_cast<size_t>(end - begin);
    if (length > MAX_LIST_LENGTH)
        error_smallListLengthExceededMax();

    List64 out = lists_getEmptyList64();

    for (const int* ptr = begin; ptr != end; ++ptr)
        out.push_back(*ptr);

    return out;
}


INLINE List64 lists_getList64(const int* elems, size_t length) {

    if (elems == nullptr && length > 0)
        error_smallListNullPtrWithPositiveLength();
    
    // no ptr necessary whgen list is empty
    if (elems == nullptr)
        return lists_getEmptyList64();

    return lists_getList64(elems, elems + length); // validates length <= MAX
}


INLINE List64 lists_getList64(std::initializer_list<int> init) {

    return lists_getList64(init.begin(), init.end());
}



/*
 * ASSERT TRIVIAL
 *
 * which doesn't really gaurantee CUDA compatibility, but may
 * catch a developer accidentally breaking compatibility
 */


static_assert(std::is_trivially_copyable_v<List64>);
static_assert(std::is_standard_layout_v<List64>);



/*
 * CONST LIST64 DECLARATION
 * 
 * Functions can accept ConstList64 (over List64) to avoid
 * a stack copy. A List64 can always be passed to a
 * function accepting a ConstList64, but a ConstList64 can never
 * be returned from a function (duh).
 */

using ConstList64 = const List64&;



#endif // LISTS_HPP
