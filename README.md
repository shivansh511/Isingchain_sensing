# Sensing with spin-s Ising Chain.

This repository contains the code files used to generate the data in the paper [arxiv:2401.14853](https://arxiv.org/pdf/2401.14853.pdf).
For the codes to function properly, please install the latest version of [Armadillo](https://arma.sourceforge.net/) C++ library for linear algebra and scietific computing and [QIClib](https://titaschanda.github.io/QIClib/documentation.html) package written in C++ for quantum computing tools based on "armadillo".
Include the path of the "QIClib" file at the top of the codes.

This repository includes two code files. Each file, in the beginning, has a global variable "spin." find it and chain it 0.5, 1.0, 1.5, 2, 2.5, and more to run the codes for different spins.


The code "qfi.cpp" calculates the quantum Fisher information of the state generated in step (I) of the protocol proposed in section (III) of the paper.
To run the code open the terminal and type "g++ qfi.cpp -o qfi.out -larmadillo && ./qfi.out 4" and press enter. Here, "4" is the length of the spin chain and can be changed to higher numbers.


File "unct.cpp" calculate the uncertainty achieved at the end of the protocol in section (III). 
To run the code open the terminal and type "g++ unct.cpp -o unct.out -larmadillo && ./unct.out 4 1 1 0.1 0". The first number "4" is the length of the spin chain, the second "1" denotes the fall-off rate and should be kept as it is for the nearest neighbor Ising chain, the third number "1" denotes the range of the interaction and should be kept to "1" for nearest neighbor chain, "0.1" is the strength of the magnetic field in the transverse direction, and final number "0" is the time from which to start calculating the and uncertainty.
Other than the length of the chain, the parameters should be left as it is for the best results.
