"""
Exercise2: Particle3D, a class to describe 3D particles

Skye Davidson
s1787851
"""
import math
import numpy as np

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    position(numpy array) - 3D position vector
    velocity(numpy array) - 3D velocity vector
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, label, x_pos, y_pos, z_pos, x_vel, y_vel, z_vel, mass):
        """
        Initialise a Particle3D instance
       
        """
       
        self.label = str(label)
        self.position = np.array([x_pos, y_pos, z_pos],float)
        self.velocity = np.array([x_vel, y_vel, z_vel],float)
        self.mass = mass
    

    def __str__(self):
        """
        Define output format.
        """
        xyz_str = "{0:8s} {1:12.8f}, {2:12.8f} {3:12.8f}".format(label, self.position[0], self.position[1], self.position[2])
        return xyz_str
        """return " x-pos =" + str(self.x_pos) + ", y-pos =" + str(self.y_pos) + ", z-pos =" + str(z_pos)

        return  "m = " + str(self.mass)"""

    
    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        return (0.5)*(self.mass)*np.dot(self.velocity,self.velocity)

    
    @staticmethod
    def kinetic_energy_list(particle_list)          """particle_list defined in main of vv"""
        ke = 0.0
        for particle in particle_list:
            ke += particle.kinetic_energy()
        return ke
        

    # Time integration methods
    def leap_vel(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m
        """
        self.velocity += dt*force/self.mass


    def leap_pos1st(self, dt):
        """
        First-order position update,
        r(t+dt) = r(t) + dt*v(t)
        """
        self.position += dt*self.velocity


    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)/m
        """
        self.position += dt*self.velocity + 0.5*dt**2*force/self.mass

    
    @staticmethod
    def leap_vel_list(particle_list, dt, forces):
        """ 
        Assume forces[i] is the force on particle_list[i]
        """

    @staticmethod
    def leap_pos1st_list(particle_list, dt):
        for i in range(particle_list):
            leap_pos1st(self, dt)

    @staticmethod
    def leap_pos2nd_list(particle_list, dt, forces):
        """ 
        Assume forces[i] is the force on particle_list[i]
        """



    @staticmethod
    def generate_particles(label,position,velocity,mass):
        """
        To create particles to be used in the code
        """
        particle = (label,position,velocity,mass)
        return particle
        
        
    @staticmethod
    def vector_sep(v1, v2):
        """
        Calculate the vector separation between two 3D particles
        """
        return (v2 - v1)
            

