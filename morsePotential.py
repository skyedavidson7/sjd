"""
CMod Ex2: symplectic Euler time integration of
a particle moving in a double well potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D


def force_dw(particle, a, b):
    force = -4 * a * particle.position ** 3 + 2 * b * particle.position
    return force


def pot_energy_dw(particle, a, b):
    potential = a * particle.position ** 4 - b * particle.position ** 2
    return potential


def u_m(r1: Particle3D, r2: Particle3D, d_e: float, r_e: float, alpha: float):
    r_12 = r2.position - r1.position
    magnitude_r_12 = np.linalg.norm(r_12)

    return d_e * ((1 - np.exp(-alpha * (magnitude_r_12 - r_e))) ** 2 - 1)


def f_1(r1: Particle3D, r2: Particle3D, d_e: float, r_e: float, alpha: float):

    r_12 = r2.position - r1.position
    magnitude_r_12 = np.linalg.norm(r_12)
    unit_r_12 = r_12 / magnitude_r_12

    q = -alpha * (magnitude_r_12 - r_e)

    return 2 * alpha * d_e * (1 - np.exp(q)) * np.exp(q) * unit_r_12


# Begin main code
def main():
    # Read name of output file from command line
    if len(sys.argv) != 3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file> <input file>")
        quit()
    else:
        outfile_name = sys.argv[1]
        input_file = sys.argv[2]

    # Open output file
    with open(outfile_name, "w") as outfile, open(input_file, "r") as infile:

        r_e = 1
        alpha = 1
        d_e = 1

        p1: Particle3D = Particle3D.create_from_file(infile)
        p2: Particle3D = Particle3D.create_from_file(infile)


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()

