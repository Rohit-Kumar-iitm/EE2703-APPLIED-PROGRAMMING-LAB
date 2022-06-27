"""
            EE2703 Applied Programming Lab - 2022
            Assignment 4: Fourier Approximations 
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 28-02-2022
"""

# Importing required modules 
from pylab import * # Importing pylab
from scipy import integrate # To integrate
from numpy import *

# Defining our functions in the code 

def Exp(x):              # Defining the function Exp(x)
    return exp(x)        # Here, returning the value of exp(x) to Exp(x), i.e., Exp(x) = exp(x)

def coscos(x):           # Defining the function coscos(x)
    return cos(cos(x))   # coscos(x) = cos(cos(x))

def u_Exp(x, k):              # Defining the function u_Exp(x, k)
    return Exp(x) * cos(k*x)  # u_Exp(x, k) = Exp(x)*cos(kx)

def v_Exp(x, k):              # Defining the function v_Exp(x, k)
    return Exp(x) * sin(k*x)  # v_Exp(x, k) = Exp(x)*sin(kx)

def u_coscos(x, k):             # Defining the function u_coscos(x, k)
    return coscos(x) * cos(k*x) # u_coscos(x, k) = coscos(x)*cos(kx)

def v_coscos(x, k):             # Defining the function v_coscos(x, k)
    return coscos(x) * sin(k*x) # v_coscos(x, k) = coscos(x)*sin(kx)

def plotExp():
    figure(1) # Naming as figure(1)
    title("Semilog plot of $e^{x}$ function") # Titling the graph 
    xlabel(r'$x\rightarrow$', size = 15) # Labelling x-axis
    ylabel(r'$e^x\rightarrow$', size = 15)	# Labelling y-axis
    semilogy(x, Exp(x), 'r', label = 'True Value')	# Plotting semilog for exp(x) function
    semilogy(x, Exp(t), '-b', label = 'Periodic Extension')	# Plotting semilog for periodic extension of exp(x) function
    semilogy(xt, Est_C_Exp, 'go', label = 'Estimated Value') # Plotting semilog for estimated exp(x) function
    grid(True) # Adding the grid lines to the plot
    legend() # Describing the elements of the graph
    show()

def plotcoscos():
    figure(2) # Naming as figure(2)		
    title("plot of cos(cos(x)) function") # Titling the graph
    xlabel(r'$x\rightarrow$',size = 15)	# Labelling x-axis
    ylabel(r'$cos(cos(x))\rightarrow$',size = 15) # Labelling y-axis
    plot(x, coscos(x),'r',label = 'True Value')	# Plotting cos(cos(x)) function
    plot(x, coscos(t),'-b',label = 'Periodic Extension') # Plotting Fourier function of cos(cos(x)) function
    plot(xt, Est_C_coscos, 'go', label = 'Estimated Value')	# Plotting estimated cos(cos(x)) function
    grid(True)  # Adding the grid lines to the plot
    legend() # Describing the elements of the graph
    show()

n = array(range(1, 52)) # Using array(range())qqqq function
def SemilogCoeffExp():
    figure(3) # Naming as figure(3)
    title("Semilog Plot of coefficients for $e^{x}$") # Titling the graph
    xlabel(r'$n\rightarrow$', size = 15) # Labelling x-axis
    ylabel('Magnitude of coefficients for $e^{x}$ ',size = 15) # Labelling y-axis	
    semilogy(n, abs(C_Exp), 'ro', label = 'Direct Integration') # semilog plot of magnitude of (Direct integration) fourier coefficients for exp(x)
    semilogy(n, abs(lstsq_C_Exp), 'go', label = 'Least Squares Approach') # semilog plot of magnitude of estimated fourier coefficients for exp(x)
    legend() # Describing the elements of the graph
    show()

def loglogCoeffExp():
    figure(4) # Naming as figure(4)
    title("Loglog Plot of coefficients of $e^{x}$")	# Titling the graph
    xlabel(r'$n\rightarrow$',size = 15) # Labelling x-axis
    ylabel('Magnitude of coefficients for $e^{x}$',size = 15) # Labelling y-axis
    loglog(n, abs(C_Exp), 'ro', label = 'Direct Integration') # loglog plot of magnitude of (Direct integration) Fourier coefficients for exp(x)
    loglog(n, abs(lstsq_C_Exp), 'go', label = 'Least Squares Approach')	# loglog plot of magnitude of estimated fourier coefficients for exp(x)
    grid(True) # Adding the grid lines to the plot
    legend() # Describing the elements of the graph
    show()

def SemilogCoeffCoscos():
    figure(5) # Naming as figure(5)
    title("Semilog Plot of coefficients for cos(cos(x))") # Titling the graph
    xlabel(r'$n\rightarrow$', size = 15) # Labelling x-axis
    ylabel('Magnitude of coefficients for cos(cos(x)) ',size = 15) # Labelling y-axis	
    semilogy(n, abs(C_coscos), 'ro', label = 'Direct Integration') # semilog plot of magnitude of (Direct integration) fourier coefficients for cos(cos(x))
    semilogy(n, abs(lstsq_C_coscos), 'go', label = 'Least Squares Approach') # semilog plot of magnitude of estimated fourier coefficients for cos(cos(x))
    grid(True) # Adding the grid lines to the plot
    legend() # Describing the elements of the graph
    show()

def loglogCoeffCoscos():
    figure(6) # Naming as figure(6)
    title("Loglog Plot of coefficients of cos(cos(x))")	# Titling the graph
    xlabel(r'$n\rightarrow$',size = 15) # Labelling x-axis
    ylabel('Magnitude of coefficients for cos(cos(x))',size = 15) # Labelling y-axis
    loglog(n, abs(C_coscos), 'ro', label = 'Direct Integration') # loglog plot of magnitude of (Direct integration) Fourier coefficients for cos(cos(x))
    loglog(n, abs(lstsq_C_coscos), 'go', label = 'Least Squares Approach')	# loglog plot of magnitude of estimated fourier coefficients for cos(cos(x))
    grid(True) # Adding the grid lines to the plot
    legend() # Describing the elements of the graph
    show()

π = np.pi # Storing the constant value(3.14159...) in π 

 # Finding the Fourier coefficients by direct integration and assigning to a matrix, say C

C_Exp = np.zeros((51, 1))	# Initialising C_Exp array to store Fourier coefficients of function exp(x)
C_Exp[0][0] = (1/(2*π))*(integrate.quad(Exp, 0, 2*π))[0] # Put x = 0, k = 0; to get a0 (Constant Fourier coefficient)		
for i in range(1, 26):									# i varies from 1 to 25
	C_Exp[(2*i)-1][0] = (1/π)*(integrate.quad(u_Exp, 0, 2*π, args=(i)))[0]  # Calculating 'an' Fourier coefficients and assigning it to odd terms of matrix(C[1], C[3]..)
	C_Exp[2*i][0] = (1/π)*(integrate.quad(v_Exp, 0, 2*π, args=(i)))[0]	# Calculating 'bn' Fourier coefficients and assigning it to even terms of matrix(C[2], C[4]..)

C_coscos = np.zeros((51,1))								# Initialising C_Exp array to store Fourier coefficients of function exp(x)
C_coscos[0][0] = (1/(2*π))*(integrate.quad(coscos, 0, 2*π))[0] # Put x = 0, k = 0; to get a0 (Constant Fourier coefficient)	
for i in range(1,26):									# i varies from 1 to 25									
	C_coscos[(2*i)-1][0] = (1/π)*(integrate.quad(u_coscos, 0, 2*π, args=(i)))[0]	# Calculating 'an' Fourier coefficients and assigning it to odd terms of matrix(C[1], C[3]..)
	C_coscos[2*i][0] = (1/π)*(integrate.quad(v_coscos, 0, 2*π, args=(i)))[0] # Calculating 'bn' Fourier coefficients and assigning it to even terms of matrix(C[2], C[4]..)


# Making up of matrix A (whcih will be used in least squares approach)

A = np.zeros((400, 51)) # Initialising A with zeros(filling 0s)
A[:, 0] = 1 # 1st column of matrix A are 1s
vector_x = linspace(0, 2*π, 401) # Defining a vector with 401 elements spaced equally between 0 and 2π
vector_x = vector_x[:-1] # Removing the last value(2π)

# Running a 'for' loop to update matrix A by filling the columns one by one
for i in range(1, 26):
    A[:, 2*i-1] = cos(i*vector_x) # Odd columns are assigned with cos(kx), x := vector_x
    A[:, 2*i] = sin(i*vector_x) # Even columns are assigned with sin(kx), x := vector_x
C_Expvector_x = Exp(vector_x)
C_coscosvector_x = coscos(vector_x)
###          ~~~GENERATING(ESTIMATING) FOURIER COEFFICIENTS THROUGH LINEAR SQUARES FITTING ~~~          ###
lstsq_C_Exp = lstsq(A, C_Expvector_x, rcond=None)[0]	# Estimating Fourier coefficients of exp(x) by lstsq function (least square approach) 
lstsq_C_coscos = lstsq(A, C_coscosvector_x, rcond=None)[0] # Estimating Fourier coefficients of cos(cos(x)) by lstsq function (least square approach)


# Calculating Transposes of matrices 
Ct_Exp = transpose(C_Exp) # Transpose of C_Exp matrix (obtained from direct integration)
Ct_coscos = transpose(C_coscos)	# Transpose of C_coscos matrix (obtained from direct integration)


###        ~~~ TOLERANCE ~~~        ###
# (i) Finding the absolute difference between the calculated and estimated Fourier coefficients
# (ii) Then finding the largest deviation in the above absolute difference

Abs_Diff_Exp = abs(lstsq_C_Exp - Ct_Exp)  # Absolute difference between calculated and estimated Fourier coefficients of exp(x)
Abs_Diff_coscos = abs(lstsq_C_coscos - Ct_coscos) # Absolute difference between calculated and estimated Fourier coefficients of cos(cos(x))

Lar_Dev_Exp = max(Abs_Diff_Exp[0]) # Largest deviation(Direct - estimated) in Fourier coefficients of exp(x) 
Lar_Dev_coscos = max(Abs_Diff_coscos[0])	# Largest deviation(Direct - estimated) in Fourier coefficients of cos(cos(x))

print("The largest deviation between coefficients for exp() is ", Lar_Dev_Exp)	#printing largest deviation value for fourier coefficients of exp(x)
print("The largest deviation between coefficients for coscos() is ", Lar_Dev_coscos)	#printing largest deviation value for fourier coefficients of cos(cos(x))	


xt = linspace(0, 2*π, 400)  # Assigning 400 values equally spaced between 0 and 2π to xt
t = tile(xt, 3) # Using tile function
x = linspace(-2*π, 4*π, 1200) # Assigning 1200 values equally spaced between -2π and 4π to x  

# Estimated functions of exp(x) and cos(cos(x))
Est_C_Exp = dot(A, lstsq_C_Exp) # Approximate exp(x)
Est_C_coscos = dot(A, lstsq_C_coscos) # Approximate cos(cos(x))

# Calling All the functions to show the graphs
plotExp()

plotcoscos()

SemilogCoeffExp()

loglogCoeffExp()

SemilogCoeffCoscos()

loglogCoeffCoscos()

    















