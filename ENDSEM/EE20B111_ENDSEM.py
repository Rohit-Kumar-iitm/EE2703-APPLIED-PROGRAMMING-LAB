"""
        EE2703 Applied Programming Lab - 2022
            END SEMESTER EXAMINATION
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 12-05-2022
"""

# Importing the required libraries
import numpy as np
from pylab import *
import math 

l = 0.5 # Quarter wavelength
c = 2.9979e8 # Speed of light
pi = math.pi
u0 = 4e-7 * pi # Permeability of free space
N = 100 # Number of sections in each half section of the antenna 
Im = 1.0 # Current injected into antenna
a = 0.01 # Radius of the wire
lamda = l*4.0 # Wavelength
f = c/lamda # Frequency
k = (2*pi)/lamda # Wave number
dz = l/N # Spacing of the current sampless

""" QUESTION - 1 """
# Calculating z vector and printing it
z = linspace(-l, l, 2*N + 1)
print(z.round(2))

# Calculating I vector(standard expression) - the actual I
I0 = array(np.zeros(2*N + 1)) # Initialisation

# Using the given expressions
I0[0:N] = Im * np.sin(k*(l + z[0:N])) # For -l < z < 0
I0[N:2*N + 1] = Im * np.sin(k*(l - z[N:2*N + 1])) # For 0 < z < l 

# Applying the given boundary condititons
I0[N] = Im
I0[0] = 0
I0[2*N] = 0 
# print(I0.round(2))

# Calculating u vector
u_index = array(range(1, 2*N))

# Removing the middlemost element
u_index = delete(u_index, N - 1)

u = z[u_index] # Final "u" vector
# print(u.round(2))

# Calculating J vector - Actual
J0 = I0[u_index] # Excluding the first, middle and extreme values of I0 using the array u
# print(J0.round(2))

""" QUESTION - 2 """
# Function to compute and return matrix M, H_phi
def Q2_make_matrices(N, J):
    M = (1/(2*pi*a))*(identity(2*N - 2, dtype = float))
    H = dot(M, J)    
    return M, H     

M, H = Q2_make_matrices(N, J0) # Getting the matrix M
# print(M.round(2))

""" QUESTION - 3 """
# Computing vectors Rz, Ru and matrices PB, P

def Q3_compute_vectors(N, r):
    zi, zj = meshgrid(z, z) # Returns coordinate matrices from coordinate vectors
    unity_matrix = ones([2*N + 1, 2*N + 1]) 
    Rz = sqrt((r**2*unity_matrix) + (zi - zj)**2) # R^2 = (r^2) + (zi - zj)^2 

    ui, uj = meshgrid(u, u)
    unity_matrix = ones([2*N - 2, 2*N - 2])
    Ru = sqrt((r**2*unity_matrix) + (ui - uj)**2) 

    return Rz, Ru

Rz, Ru = Q3_compute_vectors(N, a) # As we are evaluating at r = a
# print(Rz.round(2))
# print(Ru.round(2))

def Q3_make_matrices(N, r):
    
    RiN = Rz[N] 
    RiN = delete(RiN, [0, N, 2*N], 0) # Removing the three elements (first, middle, last) 
    PB = ((u0/(4*pi))*(exp(-1j*k*RiN))*(dz/RiN))
    
    Pij = ((u0/(4*pi))*(exp(-1j*k*Ru))*(dz/Ru)) # Using the given expressions
    return Pij, PB

Pij, PB = Q3_make_matrices(N, a) # At r = a   
# print((Pij*1e8).round(2))
# print((PB*1e8).round(2))

""" QUESTION - 4 """
# Computing the matrices Qij and QB

Qij = ((-Pij*a)/u0)*(((-1j*k)/Ru) - (1/Ru**2))

RiN = Rz[N] 
RiN = delete(RiN, [0, N, 2*N], 0)  

QB = ((-PB*a)/u0)*(((-1j*k)/RiN) - (1/RiN**2))
# print(Qij.round(2))
# print(QB.round(2))

""" QUESTION - 5 """
# Solving for the current vector, I and comparing with the actual expression
#M1 = identity(2*N-2)/(2*pi*a)
Inverse = linalg.inv(M - Qij) # Calculating the inverse of (M-Q)
J = dot(Inverse, QB) # As J = [(M-Q)^-1]QB.Im

# Finding I(expected value of current)
# Adding the three values given in question
I = zeros(2*N + 1, dtype = complex)
I[1:N] = J[0:N-1]
I[N+1:2*N] = J[N-1:2*N-1]
I[N] = Im
# print(I)

# Plotting the required graph
figure()
plot(z, I0, label = "ACTUAL", color = 'r')
plot(z, I, label = "ESTIMATED", color = 'b')
title("Graph of Antenna current", fontsize = 13)
xlabel("z → ", fontsize = 13)
ylabel("Value of current(in A) →", fontsize = 13)
legend()
grid(True)
show()





           





        


    

