 
Coded by Bernat Corominas-Murtra at the IST Austria (2020)
 
Details of the models and motivation are found in: 

  N. I. Petridou, B. Corominas-Murtra, C.-P. Heisenberg and E. Hannezo
  
  "Rigidity percolation uncovers the structural basis of embryonic tissue phase transitions"
  Accepted in Cell (2021) 
  
The author wants to thank Réka Kórei for early versions of the deformation functions
 
********Contents*********
 
*****Headers:

Lattice.py 

Constructing triangluar lattices with side length L and parameters:

    Removing a fraction of nodes

    Removing a fraction of links/obtaining a desired average connectivity

    Removing links in a correlated way --see paper for motivation and details

Deformation.py

Implements a deformation of a triangular lattice made of springs with elastic constant k 
constructed using Lattice.py
   Lateral deformation/stress applied at the bottom layer
   Vertical deformation/stress applied at the bottom layer
   Energy release due to spontaneous breaking/healing of springs --simulating viscosity, see paper

******Example

Linear_Deformation_Shear.py

It constructs a lattice using Lattice.py. Using Deformation.py, it returns 
the transmission of strain through the whole network of springs due to a small shear deformation 
applied at the top layer as a function of:
  Different --increasing-- average connectivities
  Different life times of bonds --from elastic to low viscosity
A plot of the stress sensed by the bottom layer is provided
