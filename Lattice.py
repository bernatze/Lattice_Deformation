#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:44 2019

@author: bernat
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import pebble

######################################################################
#   square lattice of length L
#
#   Returns a square lattice: 
#   if g=Regular_Lattice(L)
#    g[0]=lattice       //The topology of the graph, nodes sequentially labelled
#    g[1]=Geo_Lattice   //The links in their geometrical coorcinates 
#                         --start point, end point
#    g[2]=list_Node_Geo //Geolocation of the nodes (x,y), position in the array, 
#                         the label used in Geo_Lattice 
######################################################################

def Regular_Lattice(L):

    #list of nodes

    list_Node_Geo=[]

    for i in range (0,L):
        if (i%2)==0:
            for k in range (0,L):
                list_Node_Geo.append([i,2*k])
        else:
            for k in range(1,L):
                list_Node_Geo.append([i,2*k-1])
                
    lattice=[]
    Geo_Lattice=[]
    for i in range(0,int(np.size(list_Node_Geo)/2)-1):
        for k in range(i+1, int(np.size(list_Node_Geo)/2)):
            if (list_Node_Geo[k][0]==list_Node_Geo[i][0]):
                if(list_Node_Geo[k][1]==list_Node_Geo[i][1]+2):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo[i], list_Node_Geo[k]])
            if (list_Node_Geo[k][0]==list_Node_Geo[i][0]+1):
                if(list_Node_Geo[k][1]==list_Node_Geo[i][1]-1):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo[i], list_Node_Geo[k]])
    
            if (list_Node_Geo[k][0]==list_Node_Geo[i][0]+1):
                if(list_Node_Geo[k][1]==list_Node_Geo[i][1]+1):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo[k], list_Node_Geo[i]])
    
    return [lattice, Geo_Lattice, list_Node_Geo]

def Regular_Lattice_Symm(L):

    #list of nodes

    list_Node_Geo=[]
    list_Node_Geo_Sq=[]
    
    for i in range (0,L):
        if (i%2)==0:
            for k in range (0,L):
                list_Node_Geo_Sq.append([2*k,np.sqrt(3)*i])
                list_Node_Geo.append([2*k, i])

        else:
            for k in range(1,L):
                list_Node_Geo_Sq.append([2*k-1,np.sqrt(3)*i])
                list_Node_Geo.append([2*k-1,i])
                
    lattice=[]
    Geo_Lattice=[]
    
    for i in range(0,int(np.size(list_Node_Geo)/2)-1):
        for k in range(i+1, int(np.size(list_Node_Geo)/2)):
            if (list_Node_Geo[i][1]==list_Node_Geo[k][1]):
                if(list_Node_Geo[i][0]==list_Node_Geo[k][0]-2):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo_Sq[i], list_Node_Geo_Sq[k]])

            if (list_Node_Geo[k][1]==list_Node_Geo[i][1]+1):
                if(list_Node_Geo[k][0]==list_Node_Geo[i][0]-1):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo_Sq[i], list_Node_Geo_Sq[k]])
    
            if (list_Node_Geo[k][1]==list_Node_Geo[i][1]+1):
                if(list_Node_Geo[k][0]==list_Node_Geo[i][0]+1):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo_Sq[i], list_Node_Geo_Sq[k]])
    
    return [lattice, Geo_Lattice, list_Node_Geo_Sq]

def Regular_Lattice_Symm_Noise(L, noise):

    #list of nodes

    list_Node_Geo=[]
    list_Node_Geo_Sq=[]
    
    for i in range (0,L):

        if (i%2)==0:
            for k in range (0,L):
                P1=random.random()
                P2=random.random()
                a=noise*random.random()
                b=noise*random.random()

                if P1>0.5:
                    a=-a
                if P2>0.5:
                    b=-b
                list_Node_Geo_Sq.append([b+2*k, a+np.sqrt(3)*i])
                list_Node_Geo.append([2*k,i])

        else:
            for k in range(1,L):
                P1=random.random()
                P2=random.random()
                a=noise*random.random()
                b=noise*random.random()

                if P1>0.5:
                    a=-a
                if P2>0.5:
                    b=-b

                list_Node_Geo_Sq.append([b+2*k-1,a+np.sqrt(3)*i])
                list_Node_Geo.append([2*k-1,i])
                
    lattice=[]
    Geo_Lattice=[]
    
    for i in range(0,int(np.size(list_Node_Geo)/2)-1):
        for k in range(i+1, int(np.size(list_Node_Geo)/2)):
            if (list_Node_Geo[i][1]==list_Node_Geo[k][1]):
                if(list_Node_Geo[i][0]==list_Node_Geo[k][0]-2):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo_Sq[i], list_Node_Geo_Sq[k]])

            if (list_Node_Geo[k][1]==list_Node_Geo[i][1]+1):
                if(list_Node_Geo[k][0]==list_Node_Geo[i][0]-1):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo_Sq[i], list_Node_Geo_Sq[k]])
    
            if (list_Node_Geo[k][1]==list_Node_Geo[i][1]+1):
                if(list_Node_Geo[k][0]==list_Node_Geo[i][0]+1):
                    lattice.append([i,k])
                    Geo_Lattice.append([list_Node_Geo_Sq[i], list_Node_Geo_Sq[k]])
    
    return [lattice, Geo_Lattice, list_Node_Geo_Sq]

def Remove_Node(kill, lattice, Geo_Lattice, list_Node_Geo):

    if(len(list_Node_Geo)==0):
        print ("Empty graph!")
        return 

    else:
        del list_Node_Geo[kill]
 
        i=0 
        while i<len(lattice):            
#        for i in range (0,len(lattice)):     
            if  (i<len(lattice)):
                if ((lattice[i][0]==kill) or (lattice[i][1]==kill)):
                    del lattice[i]
                    del Geo_Lattice[i]
                else:
                    i=i+1
            else:
                break
        
        return [lattice, Geo_Lattice, list_Node_Geo]

def Remove_Random_Node(lattice, Geo_Lattice, list_Node_Geo):

    if(len(list_Node_Geo)==0):
        print ("Empty graph!")
        return

    else:
        random.seed()
        kill=[]
        kill=random.randint(0,len(list_Node_Geo)-1)
        g=Remove_Node(kill, lattice, Geo_Lattice, list_Node_Geo)
        return g

def Remove_N_Random_Nodes(N,lattice, Geo_Lattice, list_Node_Geo):

    g=[lattice, Geo_Lattice, list_Node_Geo]
    
    for i in range(0,N):
        if len(g[2])>0:
            g=Remove_Random_Node(g[0],g[1],g[2])
        else:
            print ("Empty graph!")
            g=-1
            break
        
    return g

def Remove_Link(kill, lattice, Geo_Lattice, list_Node_Geo):

    if(len(Geo_Lattice)==0):
        print ("Empty graph!")
        return 

    if(kill>len(Geo_Lattice)):
        print ("Link ID larger than number of links!")
        return 
        
    else:        
        
        del Geo_Lattice[kill]
        del lattice[kill]
        
        return [lattice, Geo_Lattice, list_Node_Geo]

def Remove_Random_Link(lattice, Geo_Lattice, list_Node_Geo):

    if(len(Geo_Lattice)==0):
        print ("Empty graph!")
        return

    else:
        random.seed()
        kill=[]
        kill=random.randint(0,len(Geo_Lattice)-1)
        g=Remove_Link(kill, lattice, Geo_Lattice, list_Node_Geo)
        return g

def Remove_N_Random_Links(N,lattice, Geo_Lattice, list_Node_Geo):

    g=[lattice, Geo_Lattice, list_Node_Geo]
    
    for i in range(0,N):
        if len(g[2])>0:
            g=Remove_Random_Link(g[0],g[1],g[2])

        else:
            print ("Empty graph!")
            g=-1
            break
    
    return g

def Extract_Nodes_Geo(Geo_lattice):
    AA=[]
    for i in range(0,len(Geo_lattice)):
        AA.append(Geo_lattice[i][0])
        AA.append(Geo_lattice[i][1])
    B=np.unique(AA, axis=0)
    return B
    

def Blastocyst_Model(L,Fluid_Fraction,Average_degree):
    g=[]
    Nodes_Pre=0
    g=Regular_Lattice_Symm_Noise(L, 0.45)
    #g=Regular_Lattice_Symm(L)
    Nodes_Pre=len(g[2])
    Nodes_Post=[]
    Frac=0
    Frac=int(Nodes_Pre*Fluid_Fraction)
    
    g=Remove_N_Random_Nodes(Frac,g[0], g[1], g[2])
    Nodes_Post=Extract_Nodes_Geo(g[1])
    Nodes=len(Nodes_Post)

    if len(g[2])>0:
        while 2*len(g[0])/Nodes>Average_degree:
            g=Remove_N_Random_Links(1,g[0], g[1], g[2])

    return g


def Blastocyst_Model_point(L,Fluid_Fraction,Average_degree):
    g=[]
    Nodes_Pre=0
    g=Regular_Lattice_Symm(L)
    Nodes_Pre=len(g[2])
    Nodes_Post=[]
    Frac=0
    Frac=int(Nodes_Pre*Fluid_Fraction)

    g=Remove_N_Random_Nodes(Frac,g[0], g[1], g[2])
    Nodes_Post=Extract_Nodes_Geo(g[1])
    Nodes=len(Nodes_Post)

    kill1=int(len(g[0])-L/2)
    kill2=int(len(g[0])-L/2)+1
    
    g=Remove_Link(kill1, g[0],g[1], g[2])
    g=Remove_Link(kill2, g[0],g[1], g[2])
        
    if len(g[2])>0:
        while 2*len(g[0])/Nodes>Average_degree:
            g=Remove_N_Random_Links(1,g[0], g[1], g[2])

    return g



def Blastocyst_Model_correlated(L,p_correlation,Average_degree):
    
    g=[]
    g=Regular_Lattice_Symm_Noise(L,0.2) ### warning, second term noise in position
    Nodes=len(g[2])
    
    mediatrices=[]
    for i in range(0, len(g[1])):
        Med=[]
        Med.append((np.array(g[1][i][0])+np.array(g[1][i][1])/2))
        mediatrices.append(Med)
        
    if len(g[2])>0:
        
        Length=len(g[0])
        Fake=g[0]
        Remain=list(range(Length))     

        while 2*Length/Nodes>Average_degree:
                        
            index= int(random.random()*len(Remain))            
            kill=Remain[index]

            Tr=1
            while Tr==1 and (2*Length/Nodes>Average_degree):
                List_links=[]
                Fake[kill][0]=-1
                Length=Length-1

                for j in range (0, len(Remain)):
                    if kill==Remain[j]:
                        del Remain[j]
                        break
                                    
                for i in range(0, len(mediatrices)):
                    
                    d1= mediatrices[kill][0][0]-mediatrices[i][0][0]
                    d2= mediatrices[kill][0][1]-mediatrices[i][0][1]
                    d=np.sqrt(d1*d1+d2*d2)
                    if (d<2) and (Fake[i][0]>-0.5): ###correlation length, citical
                        List_links.append(i) 

                P=random.random()

                ll=len(List_links)
                
                if len(List_links)>0 and (P<p_correlation):

                    while Tr==1 and (len(List_links)>0) and (2*Length/Nodes>Average_degree):
                        
                        R_link=int(random.random()*len(List_links))
                        Fake[List_links[R_link]][0]=-1
                        Length=Length-1
                        
                        for j in range (0, len(Remain)):
                            if List_links[R_link]==Remain[j]:
                                del Remain[j]
                                break
  
                        del List_links[R_link]

                        if ll>1 and len(List_links)==1:
                            Y=0
                            for w in range (0, len(Remain)):
                                if Remain[w]==List_links[0]:
                                    Y=1
                                    break

                            if P>p_correlation:
                                if Y==1:
                                    kill=List_links[0]
                                    break
                                else:
                                    Tr=0
                            else:
                                Tr=0                                
                        P=random.random()
                        if P>p_correlation:
                            Tr=0
                        if len(List_links)==0:
                            Tr=0
                else: 
                    Tr=0

    g1=[]
    g2=[]
    
    for i in range(0, len(Fake)):
        if Fake[i][0]>-1:
            g1.append(g[0][i])
            g2.append(g[1][i])

    return [g1, g2]
                
def plot_lattice(Node_Geo, GeoNet):
    plt.figure()
    plt.axis('off')
    plt.axis('equal')
    
    for i in range(0,len(GeoNet)):#Plot of the whole graph
        coord1=[]
        coord2=[]
        coord1=np.array(GeoNet[i][0])
        coord2=np.array(GeoNet[i][1])
        plt.plot([coord1[0],coord2[0]],[coord1[1],coord2[1]], '-k', alpha=.2)
    
    for i in range(0,len(Node_Geo)):#Plot of the Nodes        
        plt.plot(Node_Geo[i][0],Node_Geo[i][1], 'o', fillstyle='full', markeredgecolor='black', markerfacecolor='white')
        
    plt.show()
    
    return
    
def Max_connectivity(N_data):
    Full_connectivity=(8+16*(np.sqrt(N_data)-2)+6*(np.sqrt(N_data)-2)*(np.sqrt(N_data)-2))/N_data
    return Full_connectivity
    

def size_Giant_Rigid_Cluster(g): ###Using pebble library!!
    lattice=g[0]                
    G = pebble.lattice() # initialize a lattice
    bond=[]
                    
    for i in range(0,len(lattice)):
        if (lattice[i][0]<lattice[i][1]):            
            bond.append((lattice[i][0],lattice[i][1]))
            
        for i in bond:
            # add bond
            G.add_bond(i[0],i[1])
                     
            ############# Rigid cluster decomposition ##########
        
    G.decompose_into_cluster()
         
    size_data=0
            
    for key,value in G.cluster['index'].items():
            
        SubG=[]
        sizeSubG=0
                   
        for i in value:
            SubG.append(i)
                
        sizeSubG=np.size(np.unique(SubG))

        if sizeSubG>size_data:
            size_data=sizeSubG

    return size_data


