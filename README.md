# Impeding-waves-on-PML-
This repository contains the solution of the Maxwells equation for Electromagnetic Waves Impeding on a PML given by an spherical shell in the center of a 40 x 40 box.

The project is not a developt package. Plotting is done manually I add some scripts I use and modify to get the plots you find in the plots directory. The code only provides the E_field.txt and the H_fied.txt data as well as the amplitud evaluated at $x = 20$, $z = 20$ and $y$. Manual filter must be done for further visualization. I suggest doing it with gnuplot or with sh language of preference. 

## Compilation of the code
On Mac

`clang++ -Wall -pedantic-errors -fopenmp -O3 -o output.exe`

On Linux based systems with GCC

`g++ -Wall -pedantic-errors -fopenmp -o output.exe`

Ouputs

- On the directory data you will find the electric field E solved on space as well as the magnetic field H as well as a amplitud data of three frequencies using only the $y$ dependency of the electric field.
