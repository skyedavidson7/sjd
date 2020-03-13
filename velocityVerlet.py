"""
velocity Verlet time integration of
a particle moving in a Lennard-Jones potential.

Produces plots of the particle separation
and its energy, both as function of time. Also
saves both to file.

Skye Davidson
s1787851
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D
from MDUtilities import *


def force_lj(v1, v2):
    """
    Method to return the force on a particle
    in a Lennard-Jones potential.
    """
    r_12 = Particle3D.vector_sep(v1.position, v2.position)
    r_mod = np.linalg.norm(r_12)      

    cutoff = input("Cutoff distance is: ")

    if r_mod > cutoff:
        force_lj = 0

    else:
        force_lj = 48*((1/r_mod**14)-(1/(2*r_mod**8)))*r_12

    
    return force_lj

def forces_list(particle_list):
    """
    Method to create a list for the total forces on each particle
    """
    n = len(particle_list)
    forces_list = np.zeros([n,3])

    for i in range(n):              
        for j in range(i+1,n):
            force_ij = force_lj(particle_list[i], particle_list[j])
            forces_list[i] += force_ij
            forces_list[j] -= force_ij

    return forces_list
"""
def separation_list(particle_list, box_size, v1):
    "Method to create a list for the particle separations"
    n = len(particle_list)
    sep_vec = np.zeros([n,n,3])
    for i in range(n):
        for j in range(i+1,n):
            rij = pbc.mic(particle_list[j].position - particle_list[i].position,
                            box_size)
            sep_vec[i,j] = rij
            sep_vec[j,i] = -rij

    sep_mod = np.linalg.norm(sep_vec, axis=2)
    return sep_vec, sep_mod

    separation_list = []"""

def pot_energy_lj(v1, v2):
    """
    Method to return potential energy 
    of particle in Lennard-Jones potential
    """
    r_12 = Particle3D.vector_sep(v1.position, v2.position)
    r_mod = np.linalg.norm(r_12)                

    if r_mod > 2.5:
        potential = 0

    else:
        potential = 4*((1/r_mod**12)-(1/r_mod**6))

    return potential

def pot_energy(particle_list):
    """Method to create a list for potential energy of each particle"""
    n = len(particle_list)
    pot_energy = 0.0

    for i in range(particle_list):
        for j in range(i+1, n):
            pot_energy += pot_energy_lj(particle_list[i], particle_list[j])

    return pot_energy
        

# Begin main code
def main():
  
    # Read name of output file from command line            """do we need this"""
    if len(sys.argv)!=3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file1> <output file2> ")
        quit()
    else:
        outfile_name1 = sys.argv[1]
        outfile_name2 = sys.argv[2]

    # Open output files
    outfile1 = open(outfile_name1, "w")
    outfile2 = open(outfile_name2, "w")

    
   
    # Set up simulation parameters
    dt = 0.01
    numstep = 2000
    time = 0.0
    temp = 40
    rho = 6
    n = 4
    energy_list = []            #Initialise energy list
    energy = 0.0
    time_list = [time]          #Initialise time list
    particle_list=[]            #Initialise particle list
  
    #Create a new particle
    new_particle = Particle3D(label, position, velocity, mass)

    position = np.zeros(3)
    velocity = np.zeros(3)

    #initialize a set of P3D objects
    for i in range (n):
        new_particle = Particle3D(str(index+1), position, velocity, 1.0)
        particle_list.append(new_particle)    
         
   
    # Write out initial conditions 
    #energy = Particle3D.kinetic_energy_list(particle_list) + pot_energy(particle_list)   
    #energy_list.append(energy)
    
    #print("Energy is: " + str(energy))

    # Get initial force
    force = forces_list(particle_list)

    # Start the time integration loop
    for i in range(numstep):
        for j in range(n):
            # Update particle position
            particle_list[j].leap_pos2nd(dt, force[j])
        
            # Update force
            force_new = forces_list(particle_list)
            # Update particle velocity by averaging
            # current and new forces
            for j in range(n):
                particle_list[j].leap_vel(dt, 0.5*(force[j]+force_new[j]))
           
            # Re-define force value
            force = force_new

            # Increase time
            time += dt

            #Output file
            outfile1.write(str(n)+"\nAfile to read the particle positions\n")
            for i in particle_list:
                outfile1.write(str(i)+"\n")
            
            # Output particle information
            """energy = 0.0
            energy = Particle3D.kinetic_energy_list(particle_list) + pot_energy(particle_list)
            energy_list.append(energy)"""
        

            # Append information to data lists
            time_list.append(time)
            #pos_list.append(np.linalg.norm(Particle3D.vector_sep(v1, v2)))
            

        #print(max(energy_list) - min(energy_list))

        # Post-simulation:
        # Close output file
        outfile1.close()
        outfile2.close()

    # Plot particle trajectory to screen
    """pyplot.title('Velocity Verlet: particle separation vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Particle Separation')
    pyplot.plot(time_list, pos_list)
    pyplot.show()

    # Save plot
    pyplot.savefig('VelocityVerlet_method_particle_separation_vs_time.png')"""

    # Plot particle energy to screen
    pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()

    # Save plot
    pyplot.savefig('VelocityVerlet_method_energy_vs_time.png')


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()

