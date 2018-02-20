

We first construct a quest environment, which abstracts away all the MPI and GPU preparation.
```C
QuESTEnv env;
initQuESTEnv(&env);
```

We then create a quantum register, in this case of 3 qubits.
```C
MultiQubit qubits; 
createMultiQubit(&qubits, 3, env);
```
and set it to be in the zero state.
```C
initStateZero(&qubits);
```

We're now ready to apply some gates to our qubits, which have indices 0, 1 and 2.
When applying a gate, we pass along the quantum register.
```C
hadamard(qubits, 0);
controlledNot(qubits, 0, 1);
rotateY(qubits, 2, .1);
```

Some gates allow us to specify a general number of control qubits
```C
multiControlledPhaseGate(qubits, {0, 1, 2}, 3);
```

We can specify general single-qubit unitary operations as 2x2 matrices
```C
ComplexMatrix2 u;
u.r0c0 = (Complex) {.real=.5, .imag= .5};
u.r0c1 = (Complex) {.real=.5, .imag=-.5}; 
u.r1c0 = (Complex) {.real=.5, .imag=-.5};
u.r1c1 = (Complex) {.real=.5, .imag= .5};
unitary(qubits, 0, u);
```
or more compactly, foregoing the global phase factor,
```C
Complex a, b;
a.real = .5; a.imag =  .5;
b.real = .5; b.imag = -.5;
compactUnitary(qubits, 1, a, b);
```
or even more compactly, as a rotation around an arbitrary axis on the Bloch-sphere
```C
rotateAroundAxis(qubits, 2, 3.14/2, {1,0,0});
```

We can controlled-apply general unitaries
```C
controlledCompactUnitary(qubits, 0, 1, a, b);
```
even with multiple control qubits
```C
multiControlledUnitary(qubits, {0, 1}, 2, 2, u);
```

What has this done to the probability of the state |111> = |7>?
```C
REAL prob7 = getProbEl(qubits, 7);
```

How probable is measuring our final qubit (2) in outcome `1`?
```C
REAL prob = findProbabilityOfOutcome(qubits, 2, 1);
```

Let's measure the first qubit, collapsing it to the classical 0, 1 basis
```C
measure(qubits, 0);
```
and now measure our final qubit, while also learning of the probability of its outcome.
```
REAL prob;
int outcome = measureWithStats(qubits, 2, &prob);
```

At the conclusion of our circuit, we should free up our state-vector.
```C
destroyMultiQubit(qubits, env); 
closeQuESTEnv(env);
return 0;
```

Executing all the code above simulates the below circiut.
