"""
 CMod Ex2: Particle3D, a class to describe 3D particles
"""
import typing

import numpy as np


class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    position(float) - position of particle
    velocity(float) - velocity of particle
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, label: str, mass: float, pos: np.ndarray, vel: np.ndarray):
        """
        Initialise a Particle3D instance

        :param label: label for particle
        :param pos: position as array of float [x, y, z]
        :param vel: velocity as array of float [x, y, z]
        :param mass: mass as float
        """
        self.label = label
        self.position = pos
        self.velocity = vel
        self.mass = mass

    def __str__(self):
        """
        Define output format.
        Prints the particle in an XYZ-compatible format
        <label > <x-pos > <y-pos > <z-pos >
        """
        return self.label + " " + str(self.position[0]) + " " + str(self.position[1]) + " " + str(self.position[2])

    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        velocity = np.linalg.norm(self.velocity)
        return 0.5 * self.mass * velocity ** 2

    # Time integration methods
    def leap_velocity(self, dt: float, force: np.ndarray):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as array of float (x,y,z)
        """
        self.velocity[0] += dt * force[0] / self.mass
        self.velocity[1] += dt * force[1] / self.mass
        self.velocity[2] += dt * force[2] / self.mass

    def leap_pos1st(self, dt: float):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position[0] += dt * self.velocity[0]
        self.position[1] += dt * self.velocity[1]
        self.position[2] += dt * self.velocity[2]

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as array of float (x,y,z)
        """
        self.position[0] += dt * self.velocity[0] + 0.5 * dt ** 2 * force[0] / self.mass
        self.position[1] += dt * self.velocity[1] + 0.5 * dt ** 2 * force[1] / self.mass
        self.position[2] += dt * self.velocity[2] + 0.5 * dt ** 2 * force[2] / self.mass

    @staticmethod
    def create_from_file(fp: typing.TextIO):
        """
        Static method to return a new Particle3D instance from a file.
        File format is comma delimited in following order

        1. label
        2. mass
        3. vel_x
        4. vel_y
        5. vel_z
        6. pos_x
        7. pos_y
        8. pos_z

        :param fp: File pointer for import file
        """
        # Read a line from file and split on comma
        line = fp.readline()
        split_line = line.split(",")

        # Throw error if invalid
        if len(split_line) != 8:
            raise Exception("Invalid file input, requires 8 fields delimited by commas")

        # Pull out fields from list
        label = split_line[0]
        mass = float(split_line[1])
        vel_x = float(split_line[2])
        vel_y = float(split_line[3])
        vel_z = float(split_line[4])
        pos_x = float(split_line[5])
        pos_y = float(split_line[6])
        pos_z = float(split_line[7])
        vel = np.array([vel_x, vel_y, vel_z])
        pos = np.array([pos_x, pos_y, pos_z])

        # Return new object
        return Particle3D(label, mass, pos, vel)

    @staticmethod
    def relative_vector_separation(p1, p2):
        """
        Calculates the distance between 2 given particles

        :type p1: Particle3D
        :type p2: Particle3D
        """
        return np.linalg.norm(p1.position - p2.position)
