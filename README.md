# PDE Solver DY

This PDE solver is part of a project for the Advanced Programming course. It is concerned with the solution of the homogenous 2D heat transfer equation (LaPlace’s equation) for the steady-state case. It is based on C++ and is currently supported on Linux OS.

## Installation

### Prerequisites
The program assumes that the following tools and libraries are setup properly on the Linux OS:
- Git
- Make
- CMake
- Python (v3)
```
sudo apt-get install python3-dev
```
- Python pip
```
sudo apt install python3-pip
```
- Python numpy
```
sudo apt-get install python-numpy
```
- Matplotlib
```
pip install matplotlib
```
- Python GUI
```
sudo apt-get install python3-tk
```

### Git cloning
Clone the git repository using the following code
```
git clone https://gitlab.lrz.de/00000000014AE223/pde-solver-dy.git
```
### Building and Running
Inside the `/build` folder, run `cmake ../src` command. Then inside the `/build` folder run `make` command. Finally, open the executable file using `./pde_solver` command inside the `/build` folder.

## Sprint 1 (v1.0)
### Introduction
The solver implementation in this version is limited to a simple rectangular 2D domain in the cartesian coordinate system (i.e. _x_ and _y_). The problem definition is summarized in the figure below.

<div style="text-align:center">
    <img src="assets/images/general.png" alt="General Problem" width="450">
</div>

In this version, only Dirichlet boundary conditions are allowed for continuous homogeneous materials (only a single material for the whole domain is expected).

The user inputs and their corresponding restrictions are:
1.	Domain specification (double, double)
    -	Major dimensions: the width and height of the rectangle (_W_ and _H_, respectively). Must be greater than 0 and with compatible units.
2.	Mesh specification (int, int)
    -	Number of nodes in each direction (*n_x* and *n_y*), step size is calculated accordingly. Must be integers and greater than 1.
3.	Boundary conditions (double, double, double, double)
    -	Currently, only Dirichlet BCs are supported. Therefore, 4 temperature values are required according to the diagram above (_T1_, _T2_, _T3_, and _T4_). Temperature values are in degree Celsius.

The 5-point finite difference scheme, demonstrated in the figure below, is used to discretize the 2D heat transfer equation. The resulting linear system is solved using the LU-factorization method. 

<div style="text-align:center">
    <img src="assets/images/fd.png" alt="Finite difference scheme illustration" width="400"> <br/>
    <img src="assets/images/d2tx.png" alt="Approximation formula of second derivate in X" width="200">
    <img src="assets/images/d2ty.png" alt="Approximation formula of second derivate in Y" width="200"> 
</div>

### Outputs
After execution of the simulation with all input parameters, a text file (`results.csv`) can be found in the `/results` folder of the project containing all of the information required. The format of this __comma-separated__ file is as follows:
- Node number, 1st coordinate, 2nd coordinate, BC (1/0), T_sim, T_analytic, error

Also, if dependencies are installed correctly, a 3D plot for the simulated temperature values will be generated. Finally, the code outputs to the console the maximum, minimum, and average values of the error percentage. (if the test case option is selected).

## Test case
The simple unit test implemented for the testing of the numerical simulation accuracy is defined as special case of the general problem presented before. Only two temperature values (_T1_ and _T2_) are expected as BCs, and the inhomogeneous case is not supported. This special case is presented below:


<div style="text-align:center">
    <img src="assets/images/test.png" alt="Test case Problem" width="450">
</div>

For this problem, the analytical solution for the temperature value (_T_) for any point inside the domain is given by:

<div style="text-align:center">
    <img src="assets/images/analytical.png" alt="Analytical equation" width="400">
</div>

Results from the analytical solution is compared with the grid values obtained from the numerical solution. The error percentage for each node is calculated and the maximum, minimum, and average values are reported to the console.

## Example usage
In the `/build` folder, use the command `make` to build the executable file. Then, use the command `./pde_solver`.

# Test case:

Set Width to: 1

Set Height to: 1

Set Nodes in X direction to: 11

Set Nodes in Y direction to: 11

Set Run unit test to: 1

Set Temperature T1 to: 0

Set Temperature T2 to: 100

Check `../results/results.csv` file with the details of the results.

<div style="text-align:center">
    <img src="assets/images/3dplot.jpeg" alt="Test case visualization" width="450">
</div>

# Basic case

Set Width to: 2

Set Height to: 5

Set Nodes in X direction to: 25

Set Nodes in Y direction to: 18

Set Run unit test to: 0

Set Temperature T1 to: -20

Set Temperature T2 to: 100

Set Temperature T3 to: 50

Set Temperature T4 to: 0

Check `../results/results.csv` file with the details of the results.


<div style="text-align:center">
    <img src="assets/images/basic_case.png" alt="Basic case visualization" width="450">
</div>
