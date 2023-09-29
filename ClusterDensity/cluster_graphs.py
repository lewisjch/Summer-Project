#!/usr/bin/env python3


"""

Name: Lewis Hendry   Date: 28/08/23


Python code to read in a dump file frame by frame
and do some calculations.


Programe to calculate the density of a cluster given the identified
particles forming the cluster and plot against different values of number
of proteins.


Assumes the dump file has the following format
id type x y z ix iy iz
were x  etc. are atom coordinates, and are numbers between -L/2 and L/2
     ix etc. are image flags for taking periodic boundaries into account

Assumes that atom ids start at 1 and are continuous between 1...N

***************************************************************************
NOTE : Be careful NOT to use this code for the case of a dump line with  *
       xs ys zs as in that case the atom coordinates have been scaled    *
       to lie between 0 and 1                                            *
***************************************************************************
"""

import numpy as np
import sys
import matplotlib.pyplot as pyplot


class Atom:
    """ A Class for storing atom information """

    def __init__(self):
        """ Initialise the class """
        self.id = 0                                              # id of the atom
        self.type = 0                                            # type of the atom
        self.L = np.array([0.0,0.0,0.0],dtype=np.float64)        # size of simulation box
        self.x = np.array([0.0,0.0,0.0],dtype=np.float64)        # position of the atom
        self.image = np.array([0,0,0],dtype=np.int32)            # image flags for atoms
        self.x_unwrap = np.array([0.0,0.0,0.0],dtype=np.float64) # position of the atom - unwrapped coords
        self.unwrap_flag = False

        
    def sep(self,B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.x[0]-B.x[0])**2 + 
                        (self.x[1]-B.x[1])**2 + (self.x[2]-B.x[2])**2 )

    def minus(self,B):
        """ Subtract B.x vector from self.x vector """
        AminusB = np.array([0.0,0.0,0.0],dtype=np.float64)
        for i in range(3):
      
            AminusB[i] = self.x[i] - B.x[i]
        return AminusB

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def unwrap(self):
        """ Unwraps the coordinates for periodic box to generate x_unwrap array """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.x_unwrap[j] = self.x[j] + self.image[j]*self.L[j] # unwrap
            unwrap_flag = True



def atom_list(infile,N):
    """ Read a single frame of N atoms from a dump file 
        Expects coordinates to be in range -L/2 to L/2
        DOES NOT Unwrap coordinates for periodic box """

    atoms = [Atom() for i in range(N)]
    L = np.array([0.0,0.0,0.0],dtype=np.float64)

    # read in the 9 header lines of the dump file
    # get box size
    for i in range(9):
        line = infile.readline()
        if i==1:  ## second line of frame is timestep
            timestep = np.int32(line)
        if i==5 or i==6 or i==7:   # 6th-8th lines are box size in x,y,z dimensions
            # get the box size
            line = line.split()
            L[i-5]=np.float64(line[1]) - np.float64(line[0]);

    # now read the atoms, putting them at the correct index (index=id-1)
    for i in range(N):
        line = infile.readline()
        line = line.split()
        index = int(line[0])-1  # LAMMPS atom ids start from 1, python index starts from 0
        atoms[index].id = int(line[0])
        atoms[index].type = int(line[1])
        atoms[index].L = L
        for j in range(3):
            atoms[index].x[j] = np.float64(line[j+2])
        for j in range(3):
            atoms[index].image[j] = np.int32(line[j+5])

    return atoms, L, timestep


def lines_in_file(filename):
    """ Get the number of lines in the file """

    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def cluster_list(inf):
    """ Returns list of atoms in the cluster for the first frame """
    
    # create list of type 2 protein id's in the cluster
    cluster_ids = []
    
    # read new line every time function is run
    line = inf.readline()
    line = line.split()
    
    # note that the first number in line is the frame number
    for i in range(len(line) - 1):
        cluster_ids.append(int(line[i+1]))
        
    
    return cluster_ids
        

def properties_of_cluster(num_cluster_atoms, r):
    """ Calcuate volume and density of a sphere """
    
    volume = 4/3 * np.pi * r**3
    
    return volume, num_cluster_atoms/volume


def centre_of_mass(atoms, cluster_ids):
    """ Computes centre of mass, assuming mass of all particles are equal"""
     
    # calculate centre of mass by sum of a all postions
    com = np.zeros(3,dtype=np.float64)
    for i in range(len(cluster_ids)):
        index = cluster_ids[i] - 1
        com += atoms[index].x 
    
    return com/len(cluster_ids)


def furthest_distance(atoms, cluster_ids, com):
    """ Finds furthest particle of cluster """


    # Initialise the furthest atom variable before the loop
    furthest_atom_index = cluster_ids[0] - 1
    for i in range(len(cluster_ids)):
        index = cluster_ids[i] - 1
        
        # compare every protein with the current furthest, then overwrite the furthest
        if np.linalg.norm(atoms[index].x - com) > np.linalg.norm(atoms[furthest_atom_index].x - com):
            furthest_atom_index = index
        else:
            pass

    return np.linalg.norm(atoms[furthest_atom_index].x - com)


def new_cluster(atoms, cluster_ids, max_radius, dr, com):
    """ Creates new list of cluster atoms """
    
    # find new radius given a change in radius dr
    new_radius = max_radius - dr

    # create list of ids given updated cluster 
    new_cluster_ids = []
    
    # check every protein if they lie in the new cluster radius
    for i in range(len(cluster_ids)):
        index = cluster_ids[i] - 1
        if np.linalg.norm(atoms[index].x - com) < new_radius:
            new_cluster_ids.append(cluster_ids[i])
        else:
            pass
        
    return new_cluster_ids
    

############################################################################
#### Start of the main program

def main():
    # Read inputs from command line
    # The variable sys.argv contains whatever you typed on the command line
    # or after %run on the ipython console to launch the code.  We can use
    # it to get user inputs.
    # Here we expect one things:
    #    the name of this file
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) !=  3:
        print("You left out one or more required inputs")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} <prpr list> <no. of atoms list>" )
        sys.exit(1)
    else:
        prpr = sys.argv[1].replace('[', ' ').replace(']', ' ').replace(',', ' ').split()
        Nprot = sys.argv[2].replace('[', ' ').replace(']', ' ').replace(',', ' ').split()
        
        prpr = [float(i) for i in prpr]
        Nprot = [int(i) for i in Nprot]

    for i in Nprot:
        if i < 1000:
            print("Please enter values of Nprot >= 1000")
            sys.exit(1)
        else:
            pass            


        
    Natoms = Nprot.copy()
        
    for p in range(len(prpr)):        
        avg_cluster_radii = []
        avg_cluster_densities = []
        avg_exterior_densities = []
        
        
        for N in range(len(Nprot)):
            Natoms[N] = 1000 + Nprot[N]   # this code assumes that the first Npoly atoms are in the polymer
            
            inputdumpfile = 'dump.main_prpr'+str(prpr[p])+'_prpo1.5_Nprot'+str(Nprot[N])
            inputclusterfile = 'cluster-data_prpr'+str(prpr[p])+'_Nprot'+str(Nprot[N])

            Nlines_dump = lines_in_file(inputdumpfile)  # get length of file
            Nframes = int(Nlines_dump / (Natoms[N]+9))     # there are 9 header lines in each frame 


            # create empty lists to average data later
            cluster_radii = []
            cluster_densities = []
            exterior_cluster_densities = []
    
    
            # open cluster id file
            inf_cluster = open(inputclusterfile, 'r')
    
            # set variable since first line in file is a comment and readline() will
            # now only read data every loop
            line = inf_cluster.readline()

            # open the dump file
            inf_dump = open(inputdumpfile, 'r')  
            
    
            for f in range(Nframes):
                # create list of id's of atoms in cluster
                cluster_ids = cluster_list(inf_cluster)

                # create complete list of atoms from the dump file
                atoms, box_size, frame = atom_list(inf_dump, Natoms[N])

                # calculate centre of mass
                com = centre_of_mass(atoms, cluster_ids)
            
                # calculate maximum radius for clister to vary later
                max_radius = furthest_distance(atoms, cluster_ids, com)
        
                # create new cluster list given the approximated new radius
                new_cluster_ids = new_cluster(atoms, cluster_ids, max_radius, 4, com)
        
        
                # calculate the radius and density of the reduced cluster
                new_cluster_radius = furthest_distance(atoms, new_cluster_ids, com)
                new_cluster_volume, new_cluster_density = properties_of_cluster(len(new_cluster_ids), new_cluster_radius)
                new_density_outside = (Nprot[N]-len(new_cluster_ids)) / (box_size[0] * box_size[1] * box_size[2] - new_cluster_volume)
    
        
                #add values to list to find averages late
                cluster_radii.append(new_cluster_radius)
                cluster_densities.append(new_cluster_density)
                exterior_cluster_densities.append(new_density_outside)
                
            
            
            # close the files
            inf_cluster.close()
            inf_dump.close()
    
            
    
            # initialise variables 
            avg_cluster_radius = 0.0
            avg_cluster_density = 0.0
            avg_exterior_density = 0.0
        
            for i in range(Nframes):
                avg_cluster_radius += cluster_radii[i]/Nframes
                avg_cluster_density += cluster_densities[i]/Nframes
                avg_exterior_density += exterior_cluster_densities[i]/Nframes
    
            avg_cluster_radii.append(avg_cluster_radius)
            avg_cluster_densities.append(avg_cluster_density)
            avg_exterior_densities.append(avg_exterior_density)   
            
    
            
            
        # Plot radius against Nprot for each prpr value
        pyplot.title('Radius vs Nprot (prpr: '+str(prpr[p])+')')
        pyplot.xlabel('Nprot')
        pyplot.ylabel('Radius ($\sigma$)')
        pyplot.plot(Nprot, avg_cluster_radii)
        pyplot.show()
        
        
        # Plot density of the cluster against Nprot
        pyplot.title('Density Inside vs Nprot (prpr: '+str(prpr[p])+')')
        pyplot.xlabel('Nprot')
        pyplot.ylabel('Density ($\sigma^{-3}$)')
        pyplot.plot(Nprot, avg_cluster_densities)
        pyplot.show()
        
        # Plot density outside tthe cluster against Nprot
        pyplot.title('Density Outside vs Nprot (prpr: '+str(prpr[p])+')')
        pyplot.xlabel('Nprot')
        pyplot.ylabel('Density ($\sigma^{-3}$)')
        pyplot.plot(Nprot, avg_exterior_densities)
        pyplot.show()
        
        
        
        
        
        # Plot log radius against log Nprot for each prpr value
        pyplot.title('log(Radius) vs log(Nprot) (prpr: '+str(prpr[p])+')')
        pyplot.xlabel('log(Nprot)')
        pyplot.ylabel('log(Radius)')
        pyplot.plot(np.log(Nprot), np.log(avg_cluster_radii))
        pyplot.show()
           
        # Plot log density of the cluster against log Nprot
        pyplot.title('log(Density Inside) vs log(Nprot) (prpr: '+str(prpr[p])+')')
        pyplot.xlabel('log(Nprot)')
        pyplot.ylabel('log(Density)')
        pyplot.plot(np.log(Nprot), np.log(avg_cluster_densities))
        pyplot.show()
        
        # Plot log density outside the cluster against log Nprot
        pyplot.title('log(Density Outside) vs log(Nprot) (prpr: '+str(prpr[p])+')')
        pyplot.xlabel('log(Nprot)')
        pyplot.ylabel('log(Density)')
        pyplot.plot(np.log(Nprot), np.log(avg_exterior_densities))
        pyplot.show()



# Finished!
if __name__ == "__main__":
    main()