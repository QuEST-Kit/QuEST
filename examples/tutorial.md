

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
rotateY(qubits, 2);
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
compactUnitary(qubits, 0, a, b);
```

We can controlled-apply general unitaries
```C
controlledCompactUnitary(qubits, 0, 1, a, b);
```
even with multiple control qubits
```C
multiControlledUnitary(qubits, {0, 1}, 2, 2, u);
```

