This is result reproduction

Part 1: M-ary Phase Shift Keying

It is desired to simulate an M-ary PSK system for M = 2 (BPSK), M=4 (QPSK), M=8 (8-PSK) and M=16 (16-PSK) and compare simulation results with the theoretical results. The simulation model for M = 4 is shown below. The detector will correlate the received vector with each of the M (M=4 in this case) transmitted signal vectors (i.e. perform the dot product and selects the largest). These vectors are two-dimensional signal constellations).
The noise is generated as two statistically independent, zero mean Gaussian noise (The In-phase and Quadrature phase components) with variance 2 (note No/2 = 2, where No/2 is the two- sided PSD of the AWGN). Es (symbol energy) is normalized to 1 and Es/No is controlled by varying the Noise Variance 2. We have, Es = kEb, where Eb is the Bit Energy (so when M = 2 (BPSK), Es = Eb = 1, for M = 4, Es = 1 and Eb = 1/2 and so on). The simulation is performed for 250,000 symbols at different values of Eb/No.
The simulated results are compared with the Theoretical results. Scatter plots shows the effect of the noise variance on the constellation signal. The effects of non-idealities in the “I” and the “Q” channels at the receiver are investigated by adding combinations of amplitude offsets of 1% and 5% and phase offsets of 1deg and 5deg.

Part 2: M-ary Quadrature Amplitude Modulation

The purpose of this part is to perform Simulation of 16-QAM modulation scheme that uses rectangular single grid. The sixteen points are (1, 1), (1, 3), (3, 1) and (3, 3). The random number generator is used to generate a sequence of information bits corresponding to the 16 possible 4-bit combination b1, b2, b3, b4. These 16-combinations are mapped to the 16 signal vectors, which have coordinates of {Ami, Amq}. The detector will compute the distance between the received vector r = (ri, rq} and each of the transmitted vectors and choose the signal vector that is closest to r. The Symbol and Bit Error Probabilities vs. Eb/No (not Es) for at least 250,000 symbols (each symbol represent 4-bits) are simulated and compared with the theoretical results.

For more: wasswashafik@stu.yazd.ac.ir 
