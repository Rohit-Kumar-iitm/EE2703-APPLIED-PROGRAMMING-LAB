"""
            EE2703 Applied Programming Lab - 2022
            Assignment 5: The Resistor Problem 
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 11-03-2022
"""
# Importing required modules									
from pylab import *	# Importing pylab
import sys	# Importing sys module
import mpl_toolkits.mplot3d.axes3d as p3
from numpy import *

# Initializing the parameters with their default values
Nx, Ny = 25, 25; # Assigning Nx(size along x) and Ny(size along y) as 25
radius = 8;	# Assigning radius as 8
Niter = 1500; # Assigning Niter(no of iterations) is equal to 1500 

# Parameters given through command line
if(len(sys.argv) > 1 & len(sys.argv) < 4):
	Nx = int(sys.argv[1])
	Ny = int(sys.argv[2])
	Niter = int(sys.argv[3])

# Defining the required functions
	
def Fitting_Exp(x, A, B): # Defining the function Fitting_Exp(x, A, B)
    return A*np.exp(B*x) # Returning the value to A*np.exp(B*x) 
    
def Error_fitting(x, y): # Defining the function Error_Fitting(x, y) 
    logy = np.log(y)

    vector_x = np.zeros((len(x), 2))
    vector_x[:, 0] = x
    vector_x[:, 1] = 1
    B, logA = np.linalg.lstsq(vector_x, np.transpose(logy))[0]
    return (np.exp(logA), B)	# Returning the value to (np.exp(logA), B)

# Creating the respective matrices and initializing them.
phi = zeros((Ny, Nx)) # Initialising the matrix with dimensions Nx*Ny
x = linspace(-0.5, 0.5, Nx) # Assigning values for x from -0.5 to 0.5, Nx values
y = linspace(-0.5, 0.5, Ny) # Assigning values for x from -0.5 to 0.5, Ny values
n = array(range(Niter)) # Using array(range()) function
niter = array(range(500, 1500, 1)) # Assigning values for niter from 500 to 1500
X, Y = meshgrid(x, -y) # Creating a rectangular grid out of two given one-dimensional arrays x and -y representing the Matrix indexing to X,Y.

# The coordinates can be found for points inside the given radius.
A = (X*X) + (Y*Y)
ii = where(A <= (0.35*0.35))

# Alloting the value of Phi as 1.0 at those coordinates.
phi[ii] = 1.0

# Perform the iteration and to calculate the error in the potential.
errors = empty((Niter, 1))
for k in range(Niter):
	oldphi = phi.copy()
	phi[1:-1, 1:-1] = 0.25*(phi[1:-1, 0:-2] + phi[1:-1, 2:] + phi[0:-2, 1:-1] + phi[2:, 1:-1])

# Applying the boundary conditions.
	phi[1:-1, 0] = phi[1:-1, 1]
	phi[1:-1, -1] = phi[1:-1, -2]
	phi[0, 1:-1] = phi[1, 1:-1]
	phi[ii] = 1.0
	errors[k] = (abs(phi - oldphi)).max();

# The exponent part of the error values.
c_approx_500 = lstsq(c_[ones(Niter - 500), array(range(500, Niter))], log(errors[500:]), rcond = None) # Estimating laplace by lstsq function above 500 iterations. 
A_500, B_500 = exp(c_approx_500[0][0]), c_approx_500[0][1]
print("The values of A and B for the iterations after 500 are: ", A_500, B_500)

c_approx = lstsq(c_[ones(Niter), array(range(Niter))], log(errors), rcond = None) # Estimating laplace solutions by lstsq function for all iterations.
A, B = exp(c_approx[0][0]), c_approx[0][1]
print("The values of A and B are: ", A, B)

# The current density vectors 
Jx = zeros((Ny, Nx)) # Initialisation
Jy = zeros((Ny, Nx))
Jx[1:-1, 1:-1] = 0.5*(phi[1:-1, 0:-2] - phi[1:-1, 2:])
Jy[1:-1, 1:-1] = 0.5*(phi[2:, 1:-1] - phi[0:-2, 1:-1])


# Plotting of the initial potential contour.
figure(0)
plot(ii[0]/Nx-0.48, ii[1]/Ny-0.48, 'ro', label = "V = 1")
title("Initial Potential Contour")
xlim(-0.5, 0.5)
ylim(-0.5, 0.5)
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
grid(True)
legend()

# Plotting of error vs iteration in semilog.
figure(1)
semilogy(n,errors)
semilogy(n[::50], errors[::50], 'ro')
title("Error versus iteration")
xlabel(r'$Iteration\rightarrow$', size = 15)
ylabel(r'$Error\rightarrow$', size = 15)
grid(True)

# Plotting of error vs iteration in loglog.
figure(2)
loglog(n,errors)
loglog(n[::50], errors[::50], 'ro')
title("Error versus iteration in a loglog plot")
xlabel(r'$Iteration\rightarrow$', size = 15)
ylabel(r'$Error\rightarrow$', size = 15)
grid(True)

# Plotting the best fit of error in semilog.
fig3, ax1 = plt.subplots()
ax1.semilogy(range(Niter)[::50], errors[::50], label = 'original')
ax1.semilogy(range(Niter)[::50], Fitting_Exp(range(Niter)[::50], A, B), label = 'fit1')
ax1.semilogy(range(Niter)[::50], Fitting_Exp(range(Niter)[::50], A_500, B_500), label = 'fit2')
title("Best fit for error on semilog scale ")
xlabel(r'$Iteration\rightarrow$', size = 15)
ylabel(r'$Error\rightarrow$', size = 15)
grid(True)
legend()

# Plotting the best fit of error in loglog.
fig4, ax2 = plt.subplots()
ax2.loglog(range(Niter)[::50], errors[::50], label = 'original')
ax2.loglog(range(Niter)[::50], Fitting_Exp(range(Niter)[::50], A, B), label = 'fit1')
ax2.loglog(range(Niter)[::50], Fitting_Exp(range(Niter)[::50], A_500, B_500), label = 'fit2')
title("Best fit for error on loglog scale ")
xlabel(r'$Iteration\rightarrow$', size = 15)
ylabel(r'$Error\rightarrow$', size = 15)
grid(True)
legend()

# Plotting of the contour of potential.
figure(5)
contourf(X,Y,phi)
plot(ii[0]/Nx - 0.48, ii[1]/Ny - 0.48, 'ro', label = "V = 1")
title("Contour plot of potential")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
colorbar()
grid(True)
legend()

# Plotting the surface plots of potential.
fig1=figure(6)
ax = p3.Axes3D(fig1, auto_add_to_figure = False) 
fig1.add_axes(ax)
title("The 3-D surface plot of the potential")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
surf = ax.plot_surface(X, Y, phi, rstride = 1, cstride = 1, cmap = cm.jet)
fig1.colorbar(surf)

# plotting of the current vector plot along with the potential.
figure(7)
quiver(X, Y, Jx, Jy)
plot(ii[0]/Nx - 0.48, ii[1]/Ny - 0.48, 'ro')
title("Vector plot of the current flow")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
show()