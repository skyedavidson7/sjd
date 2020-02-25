"""
Program to provide a set of funtions to return the image of a particle from a periodic box in its original box, and return the image of the particle closest to the particle

position vector of particle = v
length of side of unit cube = l

Skye Davidson
s1787851
"""

import numpy as np

#create values for user input of vector and length
i_comp = float(input("Input i component of vector:"))
j_comp = float(input("Input j component of vector:"))
k_comp = float(input("Input k component of vector:"))
l = float(input("Input length of box:"))

#create vector
vector = [i_comp, j_comp, k_comp]

#create array
v = np.array([vector])

#print out array
print("Vector v is:", v)

def pbc(v, l):
    """
    Periodic boundary condition calculation returning image of particle in original box
    """
    
    image = np.mod(v, l)
    return image


def mic(v, l):
    """
    Minimum image convention calculation returning closest image of particle
    """

    
    closest_image = (np.mod(v + l/2, l) - l/2)
    return closest_image


#test functions from pbc.py
print("Image in box is:", pbc(v, l))

print("Closest image is:", mic(v, l))

    
