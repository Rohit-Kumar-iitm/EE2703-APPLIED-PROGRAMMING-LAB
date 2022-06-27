"""
            EE2703 Applied Programming Lab - 2022
            Assignment 3: Fitting Data To Models
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 18-02-2022
"""

#Importing required libraries to accomplish the task
import numpy as np #Numerical Python
from pylab import * 
import scipy.special as sp #Scientific Python
from scipy.linalg import lstsq #module used to return least-squares solution to the matrix equation

"""
    Read the "fitting.dat" file to extract the data from it,
    Also, making sure that this file is present in the same folder as in which the python file is present
"""

try:
    data = np.loadtxt("fitting.dat")       #Loading the text into the variable "data"
except IOError:
    print('Error.Keep the fitting.dat file in same folder') #Return a error message
    exit()

std_dev = np.logspace(-1, -3, 9) #Assigning values to (arrray)Sigma (i.e., standard deviation) through logspace.
# Here, logspace used above retuns 9 numbers between 10^(-3) to 10^(-1) which are equally spaced in log scale! (in logscale, it is from -3 to -1)  

row, col = data.shape # Finding the dimensions of the data, i.e., no. of rows and columns using "shape"(a tuple)

t = np.linspace(0, 10, row) # Finding N(=no. of rows) number of time(s) (in seconds) which are equally spaced in the interval (0, 10)
# Equally distributed(:= Number of rows(here, 101)) time intervals between 0 and 10 


# Defining our function excluding the noise
def g(t, A, B): 
    return A*sp.jv(2, t) + B*t # jv(n, t) --> Bessel Function

import matplotlib.pyplot as plt 
# Now using .plt will make the graph directly


###        ~~~ Code for the graph in Q-4 ~~~          ###

for i in range(col - 1):
    # Define our legend
    plot(data[:, 0], data[:, i + 1], label = '$\sigma$' + "=" + str(np.around(std_dev[i], 3))) 
    # Here, the standard deviation is rounded off to 3 decimal places as it were expexted from sample plots.
# Here, the values on x - axis is the 1st column of the "fitting.dat" file which is nothing but the time.
# The values on y - axis is functional values of the function g(t, A, B) at the times which is already defined

plt.legend()
# Plot the above legend on our graph

# Labelling the x and y axis
xlabel(r'$t$', size = 15)
ylabel(r'$f(t)+n$', size = 15)

# Giving a title to our graph
title(r'Q4: Data to be fitted to theory')
grid(True)

# Plotting the original graph with true value
plot(t, g(t, 1.05, -0.105), label = r"True value") 
plt.legend()

# Displaying the graph
show() 

# Code for the graph in Q-5
errorbar(t[::5], data[::5, 1], std_dev[1], fmt='ro', label = r"Error bar") 
xlabel(r"$t$", size = 15)
title(r"Q5: Data points for $\sigma$ = 0.1 along with exact function") # Title of the graph
plt.legend() 
plot(t, g(t, 1.05, -0.105), label = r"True value")
plt.legend() # Plot the above legend on y - axis
show()

# Defining our matrix M , Q 
M = np.zeros((row, 2)) # Order of matrix is (no. of rows(given)) by 2 and is initialised
Q = np.zeros((row, 1)) # Initialising the matrix Q, here Q = MP , where P = (1.05  -0.105)T    {T is transpose} 
for i in range(row):
    M[i, 0] = sp.jv(2, data[i, 0]) #Filling the first column elements of matrix M
    M[i, 1] = data[i, 0] # Filling the second column elements of matrix M
    Q[i] = M[i, 0]*1.05 + M[i, 1]*(-0.105)

###          ~~~ Function to compare two column matrices A and B ~~~         ### 
def ismatrix_equal(matA, matB):
    count = 0 # count is a temporary variable used to store the no. of elements matched
    
    for i in range(0, row): # Running through all rows
        if matA[i] == matB[i]: # If the element matches                       
            count += 1 # Increment count

    if count == row: # If all elements matched
        return True # Return true
    else: # Else, return false
        return False

###         ~~~ Defining our matrix A, B(as asked) ~~~~         ###
A = linspace(0, 2, 20) # A = 0, 0.1, ....., 2
B = linspace(-0.2, 0, 20) # B = -0.2, -0.19, ....., 0

temp = g(data[:, 0], 1.05, -0.105) # temp = g(t, A0, B0) is a vector, in assignment it is fk 


###         ~~~ Calculation of epsilon[i][j] ~~~                ###
epsilon = np.zeros((len(A), len(B))) #epsilon --> "mean squared error"

# Using dual for loops , as epsilon[i][j] = (1/101)* sigma(k = 0 to 101) {fk - g(tk, Ai, Bj)}
for i in range(len(A)):
    for j in range((len(B))):
        epsilon[i,j] = np.mean(np.square(temp - g(t, A[i], B[j]))) 


# Code for the graph in Q - 8
contplot = plt.contour(A, B, epsilon, 15) 
plot(1.05, -0.105, "ro")
annotate(r"$Exact\ location$", xy=(1.05, -0.105))
plt.clabel(contplot, inline = True)
plt.xlabel(r"$A$", size = 15)
plt.ylabel(r"$B$", size = 15)
plt.title(r"Q8: Countour plot for $\epsilon_{ij}$")
show() # Plotting the contour plots      

###        ~~~ Block of code to calculate Aerror, Berror ~~~        ### 
pred = [] # Initialising the required variables
Aerror = []
Berror = []
y_true = g(t, 1.05, -0.105) # True graph 
for i in range(col - 1):
    p, resid, rank, sig = lstsq(M, data[:, i + 1])
    aerr = np.square(p[0] - 1.05) # Auxiliary variable to hold error in A in each loop
    ber = np.square(p[1] + 0.105) # Auxiliary variable to hold error in B in each loop  
    Aerror.append(aerr) # Updating the error in A
    Berror.append(ber) # Updating the error in B


# Block used below to display graph for Q10
plot(std_dev, Aerror, "ro", linestyle = "--", linewidth = 1, label = r"$Aerr$") # Legend for Aerr in graph
plt.legend()
plot(std_dev, Berror, "go", linestyle = "--", linewidth = 1, label = r"Berr") # Legend for Berr in graph
plt.legend()
grid(True)
plt.xlabel(r"Noise standard deviation") 
plt.ylabel(r"$MS Error$", size = 15)
plt.title("$Q10:Variation\ of\  error\  with\  noise$")
show()


# Block used below to display graph for Q11
plt.loglog(std_dev, Aerror, "ro")
plt.errorbar(std_dev, Aerror, np.std(Aerror), fmt="ro", label=r"$Aerr$") # Legend for Aerr in graph
plt.legend()
plt.legend()
plt.loglog(std_dev, Berror,"go")
plt.errorbar(std_dev, Berror, np.std(Berror), fmt="go", label=r"$Berr$") # Legend for Berr in graph
plt.legend()
grid(True)
plt.ylabel(r"$MS Error$", size = 15)
plt.title(r"$Q11: Variation\ of\ error\ with\ noise$")
plt.xlabel(r"$\sigma_{n}$", size = 15)
show()


# Printing the output whether Q is equal to temp or not, where Q = MP and temp = g(t, A0, B0)
if ismatrix_equal(Q, temp):
    print("Both the matrices(Q = MP and g(t, A0, B0)) are equal.")
else: 
    print("Both the matrices(Q = MP and g(t, A0, B0)) are not equal.")            
