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
    r_12 = Particle3D.vector_sep(v1, v2)
    r_12mod = np.linalg.norm(r_12)      """get rid of these dont need"""
    r_mod = r_12mod/sigma
    r = r_12/sigma
    cutoff = input("Cutoff distance is":)

    if r_mod > cutoff:
        force = 0

    else:
        force_lj = 48*((1/r_mod**14)-(1/(2*r_mod**8)))*r

    
    return force_lj

def forces_list(i, j):
    """
    Method to create a list for the total forces on each particle
    """
        forces_list = np.zeros([n,3])

        for i in range(number_particles):               """numberparticles defined in main"""
            for j in range(i+1,number_particles):
                force_ij = force_lj([vi],[vj])
                forces_list[i] += force_ij
                forces_list[j] -= force_ij

def separation_list(particle_list, box_size v1):
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

        separation_list = []
"""
        for i in range(N):
            for j in range(i+1,N):
                separation_ij = Particle3D.vector_sep([vi],[vj])
                separation_list[i] = separation_ij
                """

def pot_energy_mp(v1, v2, sigma):
    """
    Method to return potential energy 
    of particle in Lennard-Jones potential
    """
    r_12 = Particle3D.vector_sep(v1, v2)
    r_12mod = np.linalg.norm(r_12)                  """no sigma and dont need these functions"""
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

    #Create four particles
    p1 = Particle3D.__init__             """need something for position, velocity, mass - x,y,z?"""
    p2 = Particle3D.__init__
    p3 = Particle3D.__init__
    p4 = Particle3D.__init__
    #Create vectors
    v1 = p1.position
    v2 = p2.position
    v3 = p3.position
    v4 = p4.position

   
    # Set up simulation parameters
    dt = 0.01
    numstep = 2000
    time = 0.0
    temp = 40
    number_particles = 4
  
    # Write out initial conditions
    energy = p1.kinetic_energy() + p2.kinetic_energy() + p3.kinetic_energy() + p4.kinetic_energy() + pot_energy_mp(v1, v2, sigma)   "pot energy between 4 particles?"
    outfile1.write("{0:f} {1:f}\n".format(time,np.linalg.norm(Particle3D.vector_sep(v1, v2))))
    outfile2.write("{0:f} {1:f}\n".format(time,energy))
    print(energy)

    # Get initial force
    force = force_lj(v1, v2)

    # Initialise data lists for plotting later
    time_list = [time]
    pos_list = [np.linalg.norm(Particle3D.vector_sep(v1, v2))]
    energy_list = [energy]

    # Start the time integration loop
    for i in range(numstep):
        for j in range(number_particles):
            # Update particle position
            p1.leap_pos2nd(dt, force)
            p2.leap_pos2nd(dt, -1*force)
        
            # Update force
            force_new = force_lj(v1, v2)
            # Update particle velocity by averaging
            # current and new forces
            p1.leap_vel(dt, 0.5*(force+force_new))
            p2.leap_vel(dt, 0.5*(-1*force-force_new))
            
            # Re-define force value
            force = force_new

            # Increase time
            time += dt
            
            # Output particle information
            energy = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_mp(v1, v2)
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

