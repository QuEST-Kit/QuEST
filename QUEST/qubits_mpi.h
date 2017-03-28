void initQUESTEnv(QUESTEnv *env);

double calcTotalProbability(MultiQubit multiQubit);

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta);

double findProbabilityOfZero(MultiQubit multiQubit, const int measureQubit);
