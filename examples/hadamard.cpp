#include <cmath>
#include <cstdio>
#include <complex>
#include "quest.h"

void printQureg(const Qureg QR) {
  std::printf("Qureg = \n{\n");
  for (std::size_t idx = 0; idx < QR.numAmps; ++idx) {
    std::printf("\t(%g, %gi)\n", QR.cpuAmps[idx].real(), QR.cpuAmps[idx].imag());
  }
  std::printf("}\n");
  return;
}

// Simple test program which applies a Hadamard to one qubit

int main (void)
{
  initQuESTEnv();

  reportQuESTEnv();

  Qureg qr = createQureg(1);

  reportQureg(qr);
  
  initZeroState(qr);

  const qcomp HVAL = 1.0 / std::sqrt(2);
  const CompMatr1 H = getInlineCompMatr1(
    {
      {HVAL, HVAL},
      {HVAL, -HVAL}
    }
  );

  unitary(qr, 0, H);

  printQureg(qr);

  destroyQureg(qr);

  return 0;
}