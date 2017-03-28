double calcTotalProbability(Circuit circuit);

void rotateQubit(const int rotQubit, Complex alpha, Complex beta,
                Circuit *circuit);

double findProbabilityOfZero(Circuit *circuit, const int measureQubit);
