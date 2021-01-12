7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:52:23 2020

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

import numpy as np 
import random

#############################################################
#                                                           #
#                 Class node and link                       #
#                                                           #
# Everything is coded from the node/neighbours perspective  #
# Links only for the plots                                  #
#                                                           #
#############################################################

class node(object):

    def __init__(self,idx):
        
        self.idx = idx ### index
        self.force = [0,0] #### Vector of forces
        self.vel = [0,0] #### Vector of velocity
        self.isFixed = 0 ##### if a node belogs to first row is fixed
        self.isKilled = 0 #### if a node has no connections is killed
        self.Geo_position=[0.0,0.0] ##### geomtrical position of the node
        self.isTop=0 #### if the node blongs to the topmost layer
        self.neighbours=[] #### idx of the neighbours
        self.neighbours_Geo= [] #### geometrical position of the neighbours
        self.rest_lengths=[] #### rest lengths of the links connecting neighbours
        self.lengths=[] ### actual lengths of the links connecting to neighbours
        self.Delta=[] #### stretching of the links w/ resopect the rest lengths
        self.is_percolating=1 #### 1 if there is a path from node to the bottom
                              #### 0 otherwise. Only relevant for the top layer
        self.initial_position=[0,0] ### initial position of the node
        self.is_aspired=0 ### 1 if it is part of the cells of the upper layer
                          ### that they are aspired
        
class link(object):

    def __init__(self,idx):
        
        self.idx = idx ### index
        self.Topo_link=[0,0] #### Topological link
        self.isKilled = 0 #### if is killed ==1
        self.l0 = 2.0 #### rest length by default
        self.Geo_link=[[0.0,0.0],[0.0,0.0]] #### Geometric coordinates of the link

#################################################################################
#                                                                               #
#                          *****Basic functions*****                            #    
#                                                                               #    
#################################################################################

#############################################
### returns the first index of the last row #
#############################################
        
def Who_is_Top(g,L):
    
    if L%2==0:
        Top=len(g[2])-L+1
    else:
        Top=len(g[2])-L
    return Top

###############################################
### Lattice generated from Lattice library is #
### structured in node and link classes       #
###############################################

def Struct_Lattice(g,L):
    
    Top=Who_is_Top(g,L) #### index of the Top nodes
                      
    ###### Node structure
    
    Nod=[]
    for i in range(0,len(g[2])):
        a=node(i)
        a.Geo_position=g[2][i]
        a.initial_position=g[2][i]
        if g[2][i][1]==0:
            a.isFixed=1
            
        lens=[]
        Neigb=[]
        Neigb_Geo=[]
        Real_lengths=[]
        Deltas=[]
        
        for k in range(0, len(g[0])):
            
            if g[0][k][0]==i:
                lens.append(2.0)## Rest lengths set to 2.0 by default
                Neigb.append(g[0][k][1]) ### Adding first neighbour topological
                Neigb_Geo.append(g[1][k][1]) ### Adding first neighbour geometrical
            
            if g[0][k][1]==i: #### the same but checking the other column --(i,j) or (j,i)
                lens.append(2.0)
                Neigb.append(g[0][k][0])
                Neigb_Geo.append(g[1][k][0])   
        
        for k in range(0, len(Neigb_Geo)):
            L=np.array(g[2][i])-np.array(Neigb_Geo[k])
            Real_lengths.append(np.sqrt(L[0]*L[0]+L[1]*L[1]))

        if len(Neigb)==0: ### if a node has no connections is killed
            a.isKilled=1

        Deltas=abs(np.array(Real_lengths)-np.array(lens))
        a.Delta=Deltas
        a.lengths=Real_lengths
        a.rest_lengths=lens           
        a.neighbours=Neigb   
        a.neighbours_Geo=Neigb_Geo
        if a.idx>Top-1:
            a.isTop=1
            
        Nod.append(a)
        
    ###### Link structure
    
    Link=[]
    for i in range(0,len(g[1])):
        a=link(i)
        a.Topo_link=[g[0][i][0],g[0][i][1]]
        a.Geo_link=[g[1][i][0],g[1][i][1]]
        L=[]
        L=np.array(g[1][i][0])-np.array(g[1][i][1])
        a.length=np.sqrt(L[0]*L[0]+L[1]*L[1])
        a.Delta=a.l0-a.length
        if a.Delta<0.0000001:
            a.length=2.0
            a.Delta=0.0
        Link.append(a)

    return Nod, Link

#############################################################
#### check if there is a path from node to the bottom layer #
#############################################################
    
def Path_to_Bottom_layer(Nod, node_idx): 
    T=0
    visited=[]
    for k in range(0, len(Nod[node_idx].neighbours)):
        visited.append(Nod[node_idx].neighbours[k])

    if len(visited)>0:
        
        L=1
        for i in range(0, len(Nod[node_idx].neighbours)):
            if Nod[node_idx].neighbours_Geo[i][1]==0:
                L=0
                T=1
                break
 
        k=0
        while L>0: ### Breadth forst search step
            node_idx_next=visited[k]
            for i in range(0, len(Nod[node_idx_next].neighbours)):
                if Nod[node_idx_next].neighbours_Geo[i][1]==0:
                    L=0
                    T=1
                    break
                else:
                    V=1
                    for j in range(0, len(visited)):
                        if Nod[node_idx_next].neighbours[i]==visited[j]:
                            V=0
                            break
                    if V==1:
                        visited.append(Nod[node_idx_next].neighbours[i])
                
            if k<len(visited)-1:
                k=k+1
            else:
                L=0
                break
        
    Nod[node_idx].is_percolating=T
    
###################################################################
##### Force exerted over the top layer during a time lapse dt     #
##### Differential deformation of the lattice fixed at the bottom #
###################################################################
    
def Differential_force(F0,K,damping,dt, Nod, Link):
    
    for node in Nod:
        
        x1 = node.Geo_position[0]
        y1 = node.Geo_position[1]
    
        node.force[0]=0
        node.force[1]=0
        
        position=0
        
        for neighbours in node.neighbours_Geo:
            
            x2 = neighbours[0]
            y2 = neighbours[1]
            
            l0=node.rest_lengths[position]
            position=position+1
            
            lt = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dl = lt - l0 

            sinphi = abs(y2-y1)/lt
            cosphi = abs(x2-x1)/lt
                
            SpringForceY = K*dl*sinphi
            SpringForceX = K*dl*cosphi
                
            # forces in y direction on the node
            if y1 < y2 :
                node.force[1] += SpringForceY 
                
            else :
                node.force[1] += -SpringForceY
                 
            # forces in x direction on the node    
            if x1 < x2:
                node.force[0] += SpringForceX 
                
            else :
                node.force[0] += -SpringForceX
        
        if node.isTop==1:
            node.force[1]=node.force[1]+F0
                    
        if node.isFixed == 0:
                
            node.vel[0] = node.force[0]/damping
            node.vel[1] = node.force[1]/damping
                    
            x = node.Geo_position[0] + node.vel[0]*dt
            y = node.Geo_position[1] + node.vel[1]*dt
            node.Geo_position = [x,y]
    
    for node in Nod: ##update node positions
        
        L=0
        Real_lengths=[]
        for k in range (0,len(node.neighbours)):
            node.neighbours_Geo[k]=Nod[node.neighbours[k]].Geo_position
            L=np.array(node.Geo_position)-np.array(Nod[node.neighbours[k]].Geo_position)
            Real_lengths.append(np.sqrt(L[0]*L[0]+L[1]*L[1]))

        node.lengths=Real_lengths

        node.Delta=abs(np.array(node.lengths)-np.array(node.rest_lengths))
        
    for links in Link:
            
        i=links.Topo_link[0]
        j=links.Topo_link[1]
        links.Geo_link=[Nod[i].Geo_position, Nod[j].Geo_position]
        
    return Nod, Link

###################################################################
##### Force exerted over the top layer during a time lapse dt     #
##### Differential deformation of the lattice fixed at the bottom 
##### Coordinate X of the nodes is kept fixed for the sake of stability
###################################################################
    
def Differential_force_X_Fix(F0,K,damping,dt, Nod, Link):
    
    for node in Nod:
        
        x1 = node.Geo_position[0]
        y1 = node.Geo_position[1]
    
        node.force[0]=0
        node.force[1]=0
        
        position=0
        
        for neighbours in node.neighbours_Geo:
            
            x2 = neighbours[0]
            y2 = neighbours[1]
            
            l0=node.rest_lengths[position]
            position=position+1
            
            lt = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dl = lt - l0 

            sinphi = abs(y2-y1)/lt
            cosphi = abs(x2-x1)/lt
                
            SpringForceY = K*dl*sinphi
            SpringForceX = K*dl*cosphi
                
            # forces in y direction on the node
            if y1 < y2 :
                node.force[1] += SpringForceY 
                
            else :
                node.force[1] += -SpringForceY
                 
            # forces in x direction on the node    
            if x1 < x2:
                node.force[0] += SpringForceX 
                
            else :
                node.force[0] += -SpringForceX
        
        if node.isTop==1:
            node.force[1]=node.force[1]+F0
                    
        if node.isFixed == 0:
                
            node.vel[0] = node.force[0]/damping
            node.vel[1] = node.force[1]/damping
                    
            x = node.Geo_position[0] + node.vel[0]*dt
            y = node.Geo_position[1] + node.vel[1]*dt
            node.Geo_position = [x,y]
    
    for node in Nod: ##update node positions
        
        L=0
        Real_lengths=[]
        for k in range (0,len(node.neighbours)):
            node.neighbours_Geo[k]=Nod[node.neighbours[k]].Geo_position
            L=np.array(node.Geo_position)-np.array(Nod[node.neighbours[k]].Geo_position)
            Real_lengths.append(np.sqrt(L[0]*L[0]+L[1]*L[1]))

        node.lengths=Real_lengths

        node.Delta=abs(np.array(node.lengths)-np.array(node.rest_lengths))
        
    for links in Link:
            
        i=links.Topo_link[0]
        j=links.Topo_link[1]
        links.Geo_link=[Nod[i].Geo_position, Nod[j].Geo_position]
        
    return Nod, Link

###################################################################
##### Force exerted over the top layer during a time lapse dt     #
##### Differential deformation of the lattice fixed at the bottom 
##### Coordinate Y of the nodes is kept fixed for the sake of stability
###################################################################
    
def Differential_force_Y_Fix(F0,K,damping,dt, Nod, Link):
    
    for node in Nod:
        
        x1 = node.Geo_position[0]
        y1 = node.Geo_position[1]
    
        node.force[0]=0
        node.force[1]=0
        
        position=0
        
        for neighbours in node.neighbours_Geo:
            
            x2 = neighbours[0]
            y2 = neighbours[1]
            
            l0=node.rest_lengths[position]
            position=position+1
            
            lt = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dl = lt - l0 

            sinphi = abs(y2-y1)/lt
            cosphi = abs(x2-x1)/lt
                
            SpringForceY = K*dl*sinphi
            SpringForceX = K*dl*cosphi
                
            # forces in y direction on the node
            if y1 < y2 :
                node.force[1] += SpringForceY 
                
            else :
                node.force[1] += -SpringForceY
                 
            # forces in x direction on the node    
            if x1 < x2:
                node.force[0] += SpringForceX 
                
            else :
                node.force[0] += -SpringForceX
        
        if node.isTop==1:
            node.force[1]=node.force[1]+F0
                    
        if node.isFixed == 0:
                
            node.vel[0] = node.force[0]/damping
            node.vel[1] = 0
                    
            x = node.Geo_position[0] + node.vel[0]*dt
            y = node.Geo_position[1] + node.vel[1]*dt
            node.Geo_position = [x,y]
    
    for node in Nod: ##update node positions
        
        L=0
        Real_lengths=[]
        for k in range (0,len(node.neighbours)):
            node.neighbours_Geo[k]=Nod[node.neighbours[k]].Geo_position
            L=np.array(node.Geo_position)-np.array(Nod[node.neighbours[k]].Geo_position)
            Real_lengths.append(np.sqrt(L[0]*L[0]+L[1]*L[1]))

        node.lengths=Real_lengths

        node.Delta=abs(np.array(node.lengths)-np.array(node.rest_lengths))
        
    for links in Link:
            
        i=links.Topo_link[0]
        j=links.Topo_link[1]
        links.Geo_link=[Nod[i].Geo_position, Nod[j].Geo_position]
        
    return Nod, Link

###################################################################
##### Force exerted over a region of the the top layer during a time lapse dt     #
##### Differential deformation of the lattice fixed at the bottom #
###################################################################

def Differential_force_point(F0,K,damping,dt, Nod, Link, L):
        
    for node in Nod:
        
        x1 = node.Geo_position[0]
        y1 = node.Geo_position[1]
    
        node.force[0]=0
        node.force[1]=0
        
        position=0
        
        for neighbours in node.neighbours_Geo:
            
            x2 = neighbours[0]
            y2 = neighbours[1]
            
            l0=node.rest_lengths[position]
            position=position+1
            
            lt = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dl = lt - l0 

            sinphi = abs(y2-y1)/lt
            cosphi = abs(x2-x1)/lt
                
            SpringForceY = K*dl*sinphi
            SpringForceX = K*dl*cosphi
                
            # forces in y direction on the node
            if y1 < y2 :
                node.force[1] += SpringForceY 
                
            else :
                node.force[1] += -SpringForceY
                 
            # forces in x direction on the node    
            if x1 < x2:
                node.force[0] += SpringForceX 
                
            else :
                node.force[0] += -SpringForceX
                        
        if node.Geo_position[1]>Nod[len(Nod)-1].initial_position[1]-0.1:
            if node.Geo_position[0]>Nod[len(Nod)-int(L/2)].initial_position[0]-0.1:
                if node.Geo_position[0]<Nod[len(Nod)-int(L/2)+1].initial_position[0]+0.1:
                    node.force[1]=node.force[1]+F0
                    node.is_aspired=1

        
        if node.isFixed == 0 and node.is_aspired==0:
                
            node.vel[0] = node.force[0]/damping
            node.vel[1] = node.force[1]/damping
                    
            x = node.Geo_position[0] + node.vel[0]*dt
            y = node.Geo_position[1] + node.vel[1]*dt

            node.Geo_position = [x,y]

        if node.isFixed == 0 and node.is_aspired==1:# nodes in the aspiration region keep the x coordinate

            node.vel[0] = 0
            node.vel[1] = node.force[1]/damping
        
            x = node.Geo_position[0]
            y = node.Geo_position[1] + node.vel[1]*dt

            node.Geo_position = [x,y]
    
    for node in Nod: ##update node positions
        
        LL=[]
        Real_lengths=[]
        for k in range (0,len(node.neighbours)):
            node.neighbours_Geo[k]=Nod[node.neighbours[k]].Geo_position
            LL=np.array(node.Geo_position)-np.array(Nod[node.neighbours[k]].Geo_position)
            Real_lengths.append(np.sqrt(LL[0]*LL[0]+LL[1]*LL[1]))

        node.lengths=Real_lengths

        node.Delta=abs(np.array(node.lengths)-np.array(node.rest_lengths))
        
    for links in Link:
            
        i=links.Topo_link[0]
        j=links.Topo_link[1]
        links.Geo_link=[Nod[i].Geo_position, Nod[j].Geo_position]
        
    return Nod, Link

###################################################################
##### Differential relaxation with top and bottom layer fixed     #
#####  x fixed for the whole set of nodes                         #
###################################################################

def Differential_Relaxation_X_Fix(K,damping,dt, Nod, Link):
    
    for node in Nod:
        
        x1 = node.Geo_position[0]
        y1 = node.Geo_position[1]
    
        node.force[0]=0
        node.force[1]=0
        
        position=0
        
        for neighbours in node.neighbours_Geo:
            
            x2 = neighbours[0]
            y2 = neighbours[1]
            
            l0=node.rest_lengths[position]
            position=position+1
            
            lt = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dl = lt - l0 

            sinphi = abs(y2-y1)/lt
            cosphi = abs(x2-x1)/lt
                
            SpringForceY = K*dl*sinphi
            SpringForceX = K*dl*cosphi
                
            # forces in y direction on the node
            if y1 < y2 :
                node.force[1] += SpringForceY 
                
            else :
                node.force[1] += -SpringForceY
                 
            # forces in x direction on the node    
            if x1 < x2:
                node.force[0] += SpringForceX 
                
            else :
                node.force[0] += -SpringForceX
                            
        if (node.isFixed == 0) and (node.isTop==0):
                
            node.vel[0] = node.force[0]/damping
            node.vel[1] = node.force[1]/damping
                    
            x = node.Geo_position[0] + node.vel[0]*dt
            y = node.Geo_position[1] + node.vel[1]*dt
            node.Geo_position = [x,y]
    
    for node in Nod: ##update node positions
        
        L=0
        Real_lengths=[]
        for k in range (0,len(node.neighbours)):
            node.neighbours_Geo[k]=Nod[node.neighbours[k]].Geo_position
            L=np.array(node.Geo_position)-np.array(Nod[node.neighbours[k]].Geo_position)
            Real_lengths.append(np.sqrt(L[0]*L[0]+L[1]*L[1]))

        node.lengths=Real_lengths

        node.Delta=abs(np.array(node.lengths)-np.array(node.rest_lengths))
        
    for links in Link:
            
        i=links.Topo_link[0]
        j=links.Topo_link[1]
        links.Geo_link=[Nod[i].Geo_position, Nod[j].Geo_position]
        
    return Nod, Link

###################################################################
##### Differential relaxation with top and bottom layer fixed     #
#####  x fixed for the whole set of nodes                         #
###################################################################

def Differential_Relaxation_Y_Fix(K,damping,dt, Nod, Link):
    
    for node in Nod:
        
        x1 = node.Geo_position[0]
        y1 = node.Geo_position[1]
    
        node.force[0]=0
        node.force[1]=0
        
        position=0
        
        for neighbours in node.neighbours_Geo:
            
            x2 = neighbours[0]
            y2 = neighbours[1]
            
            l0=node.rest_lengths[position]
            position=position+1
            
            lt = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dl = lt - l0 

            sinphi = abs(y2-y1)/lt
            cosphi = abs(x2-x1)/lt
                
            SpringForceY = K*dl*sinphi
            SpringForceX = K*dl*cosphi
                
            # forces in y direction on the node
            if y1 < y2 :
                node.force[1] += SpringForceY 
                
            else :
                node.force[1] += -SpringForceY
                 
            # forces in x direction on the node    
            if x1 < x2:
                node.force[0] += SpringForceX 
                
            else :
                node.force[0] += -SpringForceX
                            
        if (node.isFixed == 0) and (node.isTop==0):
                
            node.vel[0] = node.force[0]/damping
            node.vel[1] = node.force[1]/damping
                    
            x = node.Geo_position[0] + node.vel[0]*dt
            y = node.Geo_position[1] + node.vel[1]*dt
            node.Geo_position = [x,y]
    
    for node in Nod: ##update node positions
        
        L=0
        Real_lengths=[]
        for k in range (0,len(node.neighbours)):
            node.neighbours_Geo[k]=Nod[node.neighbours[k]].Geo_position
            L=np.array(node.Geo_position)-np.array(Nod[node.neighbours[k]].Geo_position)
            Real_lengths.append(np.sqrt(L[0]*L[0]+L[1]*L[1]))

        node.lengths=Real_lengths

        node.Delta=abs(np.array(node.lengths)-np.array(node.rest_lengths))
        
    for links in Link:
            
        i=links.Topo_link[0]
        j=links.Topo_link[1]
        links.Geo_link=[Nod[i].Geo_position, Nod[j].Geo_position]
        
    return Nod, Link


#################################################################################
#                                                                               #
#       *****Breaking/deforming links functions*****                            #    
#                                                                               #    
#################################################################################

################################################################
##        Deformation of rest lengths at random               ##
################################################################

def Deformation_of_rest_lengths_uniform_all(Nod, p):

    for node in range(0, len(Nod)):
        for neighbour in range (0, len(Nod[node].neighbours)):
            R=random.random()
            if R<p:
                Nod[node].rest_lengths[neighbour]=Nod[node].lengths[neighbour]
                for k in range(0, len(Nod[neighbour].neighbours)):
                    if Nod[neighbour].neighbours[k]==node:
                        Nod[neighbour].rest_lengths[k]=Nod[node].lengths[neighbour]
                        break

    return Nod
 
################################################################
##   Deformation of rest lengths as a function of stretching  ##
################################################################

def Deformation_of_rest_lengths_stretched(Nod, p, Threshold):

    for node in range(0, len(Nod)):
        if len(Nod[node].neighbours)>0:
            neighbour=random.randint(0, len(Nod[node].neighbours)-1)
            for neighbour in range (0, len(Nod[node].neighbours)):
                if Nod[node].lengths[neighbour]-Nod[node].rest_lengths[neighbour]>Threshold:
                    R=random.random()
                    if R<p:
                        #print('nyec')
                        Nod[node].rest_lengths[neighbour]=Nod[node].lengths[neighbour]
                        for k in range(0, len(Nod[neighbour].neighbours)):
                            if Nod[neighbour].neighbours[k]==node:
                                Nod[neighbour].rest_lengths[k]=Nod[node].lengths[neighbour]

    return Nod

################################################################
##   Breaking with threshold of elongation                    ##
################################################################

def Deformation_of_rest_lengths_breaking(Nod, Link, pbreak, breaking):

    for node in range(0, len(Nod)):
        if len(Nod[node].neighbours)>0:
            
            neighbour=0
            for k in range (0, len(Nod[node].neighbours)):
                if Nod[node].lengths[neighbour]-Nod[node].rest_lengths[neighbour]>breaking:
                    
                    R=random.random()
                    nod_neighbour=Nod[node].neighbours[neighbour]
                    
                    if R<pbreak:
                        print('CRAC\n')
                        del Nod[node].neighbours[neighbour]
                        del Nod[node].neighbours_Geo[neighbour]
                        del Nod[node].lengths[neighbour]
                        del Nod[node].rest_lengths[neighbour]
                        
                        kill=-1
                        for j in range (0, len(Nod[nod_neighbour].neighbours)):
                            if Nod[nod_neighbour].neighbours[j]==node:                                
                                kill=j
                                break
                        
                        del Nod[nod_neighbour].neighbours[kill]
                        del Nod[nod_neighbour].neighbours_Geo[kill]
                        del Nod[nod_neighbour].lengths[kill]
                        del Nod[nod_neighbour].rest_lengths[kill]

                        for j in range(0, len(Link)):
                            if Link[j].Topo_link[0]==node and Link[j].Topo_link[1]==nod_neighbour:
                                del Link[j]
                                break
                            if Link[j].Topo_link[1]==node and Link[j].Topo_link[0]==nod_neighbour:
                                del Link[j]
                                break    
                    else:
                        neighbour=neighbour+1
                            
                else:
                    neighbour=neighbour+1

    return Nod, Link

#################################################################################
#                                                                               #
#                   *****Final deformation functions*****                       #    
#                                                                               #    
#################################################################################

################################################################
###### Purely elastic network stretched during Time iterations #
################################################################

def constant_force_y(Nod,Link,L,F0,K,Time,dt,damping):
    
    for T in range(0,Time):

        Nod, Link=Differential_force(F0,K,damping,dt, Nod, Link)

    return Nod, Link

#################################################################
###### viscous network stretched during Time iterations         #
###### Bonds break & heal with updated rest lengths with prob p #
###### X coordinate remains fixed
#################################################################

def constant_force_y_X_Fixed(Nod,Link,L,F0,K,Time,dt,damping,p):
    
    for T in range(0,Time):

        Nod, Link=Differential_force_X_Fix(F0,K,damping,dt, Nod, Link)
        Nod=Deformation_of_rest_lengths_uniform(Nod, p)

    return Nod, Link

################################################################
###### Network with irreversible stretching and breaking 
###    Force exerted from a certain level on.
################################################################

def constant_force_y_point(Nod, Link,L,F0,K,Time,dt,damping, p, pbreak, Threshold, breaking):
    
    for T in range(0,Time):

        Nod, Link=Differential_force_point(F0,K,damping,dt, Nod, Link, L)        
        Nod=Deformation_of_rest_lengths_stretched(Nod, p, Threshold)
        Nod, Link=Deformation_of_rest_lengths_breaking(Nod, Link, pbreak, breaking)

    return Nod, Link

##################################################################
###### Normal Displacement of the Top layer and then relaxation  #
##################################################################

def Displacement_Top_Layer_X_Fix(Nod, Link,L,Displacement,K,Time,dt,damping,p):
    
    for node in Nod:
        if node.isTop==1:
            node.Geo_position[1]=node.Geo_position[1]+Displacement

    for T in range(0,Time):

        Nod, Link=Differential_Relaxation_X_Fix(K,damping,dt, Nod, Link)
        Nod=Deformation_of_rest_lengths_uniform_all(Nod, p)

    return Nod, Link

##################################################################
###### Normal Displacement of the Top layer and then relaxation  #
##################################################################

def Displacement_Top_Layer_Y_Fix(Nod, Link,L,Displacement,K,Time,dt,damping,p):
    
    for node in Nod:
        if node.isTop==1:
            node.Geo_position[0]=node.Geo_position[0]+Displacement

    for T in range(0,Time):

        Nod, Link=Differential_Relaxation_Y_Fix(K,damping,dt, Nod, Link)
        Nod=Deformation_of_rest_lengths_uniform_all(Nod, p)

    return Nod, Link


