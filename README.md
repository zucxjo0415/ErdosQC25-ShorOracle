# ErdosQC25-ShorOracle
Implements a Qiskit function `shoracle` that takes two positive integer inputs $a,N$ such that $a <N$ and outputs the circuit (an oracle for Shor's algorithm) $U = U_{a,N}$ such that


$$  U |x\rangle_1 |y\rangle_n =  |x\rangle_1 |ay \mod N\rangle_n $$
 
if $x=1$ and $y<N$, or 
 
$$ U |x\rangle_1 |y\rangle_n = |x\rangle_1 |y\rangle_n $$ 
 
otherwise. Here $n =2 \lceil \log_2 N \rceil$. The function also uses 3 additional ancillas.

This function is based off Beauregard's controlled modular multiplication circuit [1], which implements modular multiplication as an appropriate sequence of additions / subtractions. The addition / subtraction operations use Draper's QFT-based adder circuit. 

The function (including all the necessary auxiliary functions) can be found in the Python script file `Shoracle.py`.

The Jupyter notebook `Shoracle.ipynb` has documented code, including a brief explanation of the construction, basic complexity analysis and demonstrations.

**References**

[1] Stéphane Beauregard, "Circuit for Shor's algorithm using 2n+3 qubits", *Quant. Inf. Comput.* 3 (2003) 2, 175-185. DOI: 10.26421/QIC3.2-8. Preprint:[arXiv:quant-ph/0205095](https://arxiv.org/pdf/quant-ph/0205095)
