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


def force_lj(v1, v2):
    """
    Method to return the force on a particle
    in a Lennard-Jones potential.
    """
    r_12 = Particle3D.vector_sep(v1, v2, sigma)
    r_12mod = np.linalg.norm(r_12)
    r_mod = r_12mod/sigma
    r = r_12/sigma

    if r_mod > 2.5:
        force = 0

    else:
        force = 48*((1/r_mod**14)-(1/(2*r_mod**8)))*r

    
    return force


def pot_energy_mp(v1, v2, sigma):
    """
    Method to return potential energy 
    of particle in Lennard-Jones potential
    """
    r_12 = Particle3D.vector_sep(v1, v2)
    r_12mod = np.linalg.norm(r_12)
    r_mod = r_12mod/sigma
    r = r_12/sigma


    if r_mod > 2.5:
        potential = 0

    else:
        potential = 4*((1/r_mod**12)-(1/r_mod**6))

    return potential


# Begin main code
def main():
  
    # Read name of output file from command line
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

    # Create parameters using data from file
    file_handle = open(input("Use Oxygen.txt or Nitrogen.txt?"), "r")

    #Create two particles
    p1 = Particle3D.file_input(file_handle)
    p2 = Particle3D.file_input(file_handle)
    
    #Create vectors
    v1 = p1.position
    v2 = p2.position

    #Defining value of sigma
    input_data = file_handle.readline()
    token = input_data.split(" ")
    sigma = float(token[0])
   
    # Set up simulation parameters
    dt = 0.01
    numstep = 2000
    time = 0.0
  
    # Write out initial conditions
    energy = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_mp(v1, v2, sigma)
    outfile1.write("{0:f} {1:f}\n".format(time,np.linalg.norm(Particle3D.vector_sep(v1, v2))))
    outfile2.write("{0:f} {1:f}\n".format(time,energy))
    print(energy)

    # Get initial force
    force = force_mp(v1, v2, sigma)

    # Initialise data lists for plotting later
    time_list = [time]
    pos_list = [np.linalg.norm(Particle3D.vector_sep(v1, v2))]
    energy_list = [energy]

    # Start the time integration loop
    for i in range(numstep):
        # Update particle position
        p1.leap_pos2nd(dt, force)
        p2.leap_pos2nd(dt, -1*force)
        
        # Update force
        force_new = force_lj(v1, v2, sigma)
        # Update particle velocity by averaging
        # current and new forces
        p1.leap_velocity(dt, 0.5*(force+force_new))
        p2.leap_velocity(dt, 0.5*(-1*force-force_new))
        
        # Re-define force value
        force = force_new

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_mp(v1, v2, sigma)
        outfile1.write("{0:f} {1:f}\n".format(time,np.linalg.norm(Particle3D.vector_sep(v1, v2))))
        outfile2.write("{0:f} {1:f}\n".format(time,energy))


        # Append information to data lists
        time_list.append(time)
        pos_list.append(np.linalg.norm(Particle3D.vector_sep(v1, v2)))
        energy_list.append(energy)

    print(max(energy_list) - min(energy_list))

    # Post-simulation:
    # Close output file
    outfile1.close()
    outfile2.close()

    # Plot particle trajectory to screen
    pyplot.title('Velocity Verlet: particle separation vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Particle Separation')
    pyplot.plot(time_list, pos_list)
    pyplot.show()

    # Save plot
    pyplot.savefig('VelocityVerlet_method_particle_separation_vs_time.png')

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
