double calcTotalProbability(Circuit circuit);

void rotateQubit(const int rotQubit,
                double aRe, double aIm, double bRe,  double bIm,
                Circuit *circuit);

double findProbabilityOfZero(Circuit *circuit, const int measureQubit);
