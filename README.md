# Sensing with spin-s Ising Chain.

For the codes to function properly, please install the latest version of "QIClib" package written in C++ for quantum computing tools.
And include the path of the "QIClib" file at the top of the codes.

This repository includes two code files. 
The code "qfi.cpp" calculates the quantum Fisher information of the state generated in step (I) of the protocol proposed in section (III) of the paper [arxiv:2401.14853](https://arxiv.org/pdf/2401.14853.pdf).
To run the code open the terminal and type "g++ qfi.cpp -o qfi.out -larmadillo && ./qfi.out 4" and press enter. Here, "4" is the length of the spin chain and can be changed to higher numbers.
