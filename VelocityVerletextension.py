
import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D
from MDUtilities import set_initial_positions
from MDUtilities import set_initial_velocities


def main():

    number_particles=4
    numstep=5
    particles_list=[]
    '''initialize a set of P3D objects'''
    for i in range (number_particles):
        particles_list[i]=Particle3D.generate_particle     """generate_particle function in P3D"""
      

    """get & store initial separations & initial forces"""

'''use mdutilities to initialize positions and velocities'''
    temp=40                 """set temperature"""
    rho=6                      "set rho"
    MDUtilities.set_initial_positions(rho,particles_list)
    MDUtilities.set_initial_velocities(temp,particles_list)

    """need box size array and store it"""

"""loop that runs at each timestep that will update the positions, velocities and forces"""
    
    """leap position"""
    """update position"""
    "update separation"
    "update forces"
    "update velocity"
    """and then do it all over again"""
    """do updates for each timestep"""

    for i in range(numstep):
        for j in range(number_particles):


"""edit so that separations and positions are not the same"""
"""change velocity verlet so that 4 particles instead of 2"""
"""every line check that doesnt give errors, don't have output file yet"""
