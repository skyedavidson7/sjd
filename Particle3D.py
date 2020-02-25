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
        self.position = np.array([x_pos, y_pos, z_pos])
        self.velocity = np.array([x_vel, y_vel, z_vel])
        self.mass = mass
    

    def __str__(self):
        """
        Define output format.
        """
        return "label =" + label
        return " x-pos =" + str(self.x_pos) + ", y-pos =" + str(self.y_pos) + ", z-pos =" + str(z_pos)

        return  "m = " + str(self.mass)

    
    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        return (0.5)*(self.mass)*((np.linalg.norm(self.velocity))**2)
        

    # Time integration methods
    def leap_velocity(self, dt, force):
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
    def file_input(file_handle):
        """
        To assign appropriate properties after reading in from a file
        """
        #read in data as list of strings and then close file
        input_data = file_handle.readline()
        token = input_data.split(" ")
        label = str(token[0])
        x_pos = float(token[1])
        y_pos = float(token[2])
        z_pos = float(token[3])
        x_vel = float(token[4])
        y_vel = float(token[5])
        z_vel = float(token[6])
        mass = float(token[7])
        
        return Particle3D(label, x_pos, y_pos, z_pos, x_vel, y_vel, z_vel, mass)
       
        
    @staticmethod
    def vector_sep(v1, v2):
        """
        Calculate the vector separation between two 3D particles
        """
        return (v2 - v1)
            

