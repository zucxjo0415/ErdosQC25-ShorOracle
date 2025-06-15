from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.circuit.library import QFT

import numpy as np

def shoracle(a,N):
    if a >= N: raise Exception("Invalid inputs!")

    # classical precomputation
    n = int(np.ceil(np.log2(N)))
    apow = [[(a*(1<<i)%N)*np.pi/(1<<j) for j in range(n+1)] for i in range(n+1)]
    Npow = [N*np.pi/(1<<j) for j in range(n+1)]
    # find ainv using Euclid's algorithm
    ainv = pow(a, -1, N)
    aipow = [[(ainv*(1<<i)%N)*np.pi/(1<<j) for j in range(n+1)] for i in range(n+1)]
    
    ctrl = QuantumRegister(size=1, name="x")
    reg0 = QuantumRegister(size=n, name="y0")
    reg1 = QuantumRegister(size=n, name="y1")
    anci = QuantumRegister(size=3, name="z")
    circ = QuantumCircuit(ctrl, reg0, reg1, anci, name="U_a")

    # entire y register + overflow ancilla
    yreg = [1+i for i in range(2*n+1)]
    # subtract N
    circ.compose(QFT(2*n+1), qubits=yreg, inplace=True)
    circ.p(-N*np.pi, anci[0])
    for k in range(1,n+1): 
        circ.p(-N*np.pi/(1<<k),reg1[n-k])        
        circ.p(-N*np.pi/(1<<(n+k)),reg0[n-k])
    circ.compose(QFT(2*n+1).inverse(), qubits=yreg, inplace=True)
    # mark last ancilla 1 iff underflow (i.e. iff y<N)
    circ.cx(anci[0],anci[2])    
    # add N back
    circ.compose(QFT(2*n+1), qubits=yreg, inplace=True)
    circ.p(N*np.pi, anci[0])
    for k in range(1,n+1): 
        circ.p(N*np.pi/(1<<k),reg1[n-k])        
        circ.p(N*np.pi/(1<<(n+k)),reg0[n-k])
    circ.compose(QFT(2*n+1).inverse(), qubits=yreg, inplace=True)

    breg = [n+1+i for i in range(n+1)]
    
    # triple-controlled +*a mod N
    circ.compose(QFT(n+1), qubits=breg, inplace=True)
    for i in range(n):
        # add 2^i*a
        circ.mcp(apow[i][0], [anci[2],ctrl,reg0[i]], anci[0])
        for k in range(1,n+1): circ.mcp(apow[i][k], [anci[2],ctrl,reg0[i]], reg1[n-k])
        
        # subtract N
        circ.cp(-Npow[0], anci[2],anci[0])
        for k in range(1,n+1): circ.cp(-Npow[k],anci[2],reg1[n-k])        
        
        # undo subtract N if not needed
        circ.compose(QFT(n+1).inverse(), qubits=breg, inplace=True)
        circ.cx(anci[0],anci[1])
        circ.compose(QFT(n+1), qubits=breg, inplace=True)
        circ.mcp(Npow[0],[anci[2],anci[1]],anci[0])
        for k in range(1,n+1): circ.mcp(Npow[k],[anci[2],anci[1]],reg1[n-k])
        
        # clean ancilla
        circ.mcp(-apow[i][0], [anci[2],ctrl,reg0[i]], anci[0])
        for k in range(1,n+1): circ.mcp(-apow[i][k], [anci[2],ctrl,reg0[i]], reg1[n-k])
        circ.compose(QFT(n+1).inverse(), qubits=breg, inplace=True)
        circ.cx(anci[2],anci[0])
        circ.cx(anci[0],anci[1])
        circ.cx(anci[2],anci[0])
        circ.compose(QFT(n+1), qubits=breg, inplace=True)
        circ.mcp(apow[i][0], [anci[2],ctrl,reg0[i]], anci[0])
        for k in range(1,n+1): circ.mcp(apow[i][k], [anci[2],ctrl,reg0[i]], reg1[n-k])
    circ.compose(QFT(n+1).inverse(), qubits=breg, inplace=True)

    for i in range(n):
        #circ.ccswap([anci[2],ctrl],reg0[i],reg1[i])
        circ.cx(reg1[i],reg0[i])
        circ.mcx([anci[2],ctrl,reg0[i]],reg1[i])
        circ.cx(reg1[i],reg0[i])
    
    # controlled inverse-*a^{-1} mod N
    circ.compose(QFT(n+1), qubits=breg, inplace=True)
    for j in range(n):
        i = n-1-j
        # inverse (clean ancilla)
        circ.mcp(-aipow[i][0], [anci[2],ctrl,reg0[i]], anci[0])
        for k in range(1,n+1): circ.mcp(-aipow[i][k], [anci[2],ctrl,reg0[i]], reg1[n-k])
        circ.compose(QFT(n+1).inverse(), qubits=breg, inplace=True)
        circ.cx(anci[2],anci[0])
        circ.cx(anci[0],anci[1])
        circ.cx(anci[2],anci[0])
        circ.compose(QFT(n+1), qubits=breg, inplace=True)
        circ.mcp(aipow[i][0], [anci[2],ctrl,reg0[i]], anci[0])
        for k in range(1,n+1): circ.mcp(aipow[i][k], [anci[2],ctrl,reg0[i]], reg1[n-k])

        # inverse (undo subtract N if not needed)
        circ.mcp(-Npow[0],[anci[2],anci[1]],anci[0])
        for k in range(1,n+1): circ.mcp(-Npow[k],[anci[2],anci[1]],reg1[n-k])
        circ.compose(QFT(n+1).inverse(), qubits=breg, inplace=True)
        circ.cx(anci[0],anci[1])
        circ.compose(QFT(n+1), qubits=breg, inplace=True)
        
        # inverse (subtract N)
        circ.cp(Npow[0], anci[2], anci[0])
        for k in range(1,n+1): circ.cp(Npow[k],anci[2],reg1[n-k])        

        # inverse (double-controlled) add 2^i*a
        circ.mcp(-aipow[i][0], [anci[2],ctrl,reg0[i]], anci[0])
        for k in range(1,n+1): circ.mcp(-aipow[i][k], [anci[2],ctrl,reg0[i]], reg1[n-k])
    circ.compose(QFT(n+1).inverse(), qubits=breg, inplace=True)
    
    return circ