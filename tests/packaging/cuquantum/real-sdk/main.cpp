#include <custatevec.h>
#include <iostream>

int main() {
    const auto runtime = custatevecGetVersion();
    std::cout << "cuStateVec header=" << CUSTATEVEC_VERSION
              << " runtime=" << runtime << '\n';
    return runtime / 10000 == CUSTATEVEC_VER_MAJOR ? 0 : 1;
}
