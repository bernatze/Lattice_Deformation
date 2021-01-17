 
# Contents

Python scripts coded by Bernat Corominas-Murtra at the IST Austria (2020).
The author wants to thank Réka Kórei for early versions of the deformation functions.
 
## Micropipette aspiration analysis

Standard ImageJ script for analysis of the micropipette aspiration time lapse movies to acquire deformation rates.

    PipetteAssayMeasurement.txt

## Python scripts

### Header 1:

    Lattice.py 

Constructing triangluar lattices with side length L and parameters:

    Removing a fraction of nodes

    Removing a fraction of links/obtaining a desired average connectivity

    Removing links in a correlated way --see paper for motivation and details

### Header 2:

    Deformation.py

Implements a deformation of a triangular lattice made of springs with elastic constant k 
constructed using Lattice.py:
 
    Lateral deformation/stress applied at the bottom layer

    Vertical deformation/stress applied at the bottom layer

    Energy release due to spontaneous breaking/healing of springs --simulating viscosity, see paper

### Example

    Linear_Deformation_Shear.py

It constructs a lattice using Lattice.py. Using Deformation.py, it returns 
the transmission of strain through the whole network of springs due to a small shear deformation 
applied at the top layer as a function of:

    Different --increasing-- average connectivities
   
    Different life times of bonds --from elastic to low viscosity

    A plot of the stress sensed by the bottom layer is provided
