
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 14:20:43 2020

@author: bernat
"""

#################################################################
#                                                               #
#    Coded by Bernat Corominas-Murtra at the IST Austria        #
#    (2020)                                                     #
#    N. I. Petridou, B. Corominas-Murtra, C.-P. Heisenberg      #
#    and E. Hannezo                                             #
#    "Rigidity percolation uncovers the structural basis of     #
#    embryonic tissue phase transitions"                        #
#    Accepted in Cell (2021)                                    #
#    The author wants to thank Réka Kórei for early versions    #
#    of the deformation functions                               #
#                                                               #
#################################################################

#################################################################
#                                                               #
#   Plot of the linear response of the bottom layer of          #
#   a triangular network of springs with stochastic bond        # 
#   breaking/healing with probability p x time step after       #
#   shear stress F0 applied at the topmost layer                #                                 
#                                                               #
#################################################################

import matplotlib.pyplot as plt
import numpy as np
import Lattice
import Deformation as DF

################ parameters
##### damping > dt !!!! ########
L=10  ## Number of cells at the bottom of the square lattice 20
F0=0.05 ## Force applied to the top nodes
damping=0.00003 ## Damping  0.00001
K=1 ## Elastic constant of the springs 2
dt=0.00001 ## Time differential 0.000001
Time=400 ## Time differential per frame 100
p=0. ## probability of bond healing and breaking per differential ti
iterations=3 ## Number of different probabilities of bond breaking/healing explored
replicas=1 ## Replicas per connectivity
Displacement=0.01 # 1
rge=20 ## resolution of the intervals of connectivity
#
############ END parameters

F_Global=[]
last=[]
FF_all_Nets=[]

for r in range (0, iterations):
    q=0.35/Time
    p=q*r #### p=0 elastic response! p>0 viscous response!
    print('iteration %d of %d'%(r,iterations))
    conn=[]
    Force_Top=[]
    Force_Bottom=0
    FF=[]

    for u in range (0, rge):
        c=2+3.4*u/rge #### range of connectivities explored
        conn.append(c)
        print(str(100*u/rge)+'%')

        for z in range(0,replicas): 
            g=[]  ## lattice at time=0
            g=Lattice.Blastocyst_Model(L,0.0,c) ## Generating the lattice using Lattice.py
            FF_all_Nets.append(g[0])
        
            Nod=[] ## adapting the lattice to the class DF.node to work on deformations
            Link=[]
            Nod, Link=DF.Struct_Lattice(g,L)
            
            ##### Performing the displacement of the top layer towards the X axis
            Nod, Link=DF.Displacement_Top_Layer_Y_Fix(Nod, Link,L,Displacement,K,Time,dt,damping,p)
        
            Tr=1
            k=0
            F_T=[]
            while Tr==1:
                if  Nod[k].Geo_position[1]==0: #### stress transmitted to the bottom nodes
                    F_T.append(Nod[k].force[0]) #### stress sensed by the node                  
                    k=k+1
                else:
                    Tr=0
            Force_Bottom=0
            Force_Bottom=np.mean(np.array(F_T)) ### force sensed by the bottom layer

            if z==0:
                FF.append(Force_Bottom)
            else:
                FF[u]=FF[u]+Force_Bottom

    F_Global.append(np.array(FF)/replicas)
    last.append(FF[u]/replicas)

############# Plot figure with force propagation against normalized connectivity

plt.figure(figsize=(6,6))
conn=np.array(conn)/(Lattice.Max_connectivity(L*(L-1)))

for i in range(0, iterations):
    plt.plot(conn, F_Global[i], '^',  color='orange', markersize=8, alpha=1-i/iterations) 
    plt.ylabel('Force propagation $\propto \eta$', fontsize=16)
    plt.xlabel('Normalized connectivity', fontsize=16)
plt.show()

