# Numerov1D
By using Numerov method, the program solves the Schroedinger equation for a particle moving in one-dimensional potential of the following form:
V(x) = a x^2 - 1/((x-d)^2 + b^2)^2 - 1/((x+d)^2 + b^2)^2. Energies and corresponding wave functions of the system are being found and can be saved.
The potential and boundary conditions can be easily modified. A simple parallelization with OpenMP is added in order to perform calculations 
simultaneously for different values of a parameter d.
