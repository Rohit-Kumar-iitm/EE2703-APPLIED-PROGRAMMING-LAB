"""
            EE2703 Applied Programming Lab - 2022
            Assignment 6: The Laplace Transform 
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 27-03-2022
"""

# Importing the required modulues
import scipy.signal as sp
from pylab import *


# Q-1: The time response of the spring with decay 0.5 seconds and it's plot.

Num1 = poly1d([1, 0.5])                    # X(s) = {(s + 0.5)/[(s^2 + 2.25)(s^2 + s + 2.5)]}
Den1 = polymul([1, 1, 2.5], [1, 0, 2.25]) 
Xs_1 = sp.lti(Num1, Den1)
t1, xt_1 = sp.impulse(Xs_1, None, linspace(0, 50, 1500)) 

figure(0) # Plot of time response of the system
title("Q1: x(t) for decay of 0.5 seconds") 
xlabel("t → ", fontsize = 13)
ylabel("x(t) → ", fontsize = 13)
grid(True) 
plot(t1, xt_1)
show()

# Q-2: The time response of the spring with decay 0.05 seconds and it's plot.

Num2 = poly1d([1, 0.05])                  
Den2 = polymul([1, 0.1, 2.2525], [1, 0, 2.25]) 
Xs_2 = sp.lti(Num2, Den2)
t2, xt_2 = sp.impulse(Xs_2, None, linspace(0, 50, 1500)) 

figure(1) # Plot of time response of the system
title("Q2: x(t) for decay of 0.05 seconds") 
xlabel("t → ", fontsize = 13)
ylabel("x(t) → ", fontsize = 13)
grid(True) 
plot(t2, xt_2)
show() 

# Q-3: To obtain system transfer function and also to find out the responses of output by varying ω from 1.4 to 1.6 in steps of 0.05 

H = sp.lti([1], [1, 0, 2.25])
for ω in arange(1.4, 1.6, 0.05):
	t = linspace(0, 50, 1500)
	f = cos(ω * t) * exp(-0.05 * t)
	t3, x, svec = sp.lsim(H, f, t)


	figure(2)
	plot(t3, x, label = 'ω = ' + str(ω))
	title("Q3: x(t) for different frequencies(ω range from 1.4 to 1.6)")
	xlabel("t → ", fontsize = 13)
	ylabel("x(t) → ", fontsize = 13)
	legend(loc = 'upper left')
	grid(True)
show()	

# Q-4: To solve coupled differential equations using laplace transforms and to find the time domain resonse of the functions 
# by using sp.impulse function to find time response

H4_x = sp.lti(poly1d([1, 0, 2]), poly1d([1, 0, 3, 0]))
t4, x4 = sp.impulse(H4_x, None, linspace(0, 20, 1500))
H4_y = sp.lti(poly1d([2]), poly1d([1, 0, 3, 0]))
t4, y4 = sp.impulse(H4_y, None, linspace(0, 20, 1500))  

# The plots of time responses X(t) and Y(t).

figure(3)
plot(t4, x4, label = 'x(t)') 
plot(t4, y4, label = 'y(t)')
title("Q4: x(t) and y(t) (time responses)")
xlabel("t → ", fontsize = 13)
ylabel("Functions X(t) and Y(t) → ", fontsize = 13)
legend(loc = 'upper right')
grid(True)
show()

# Q-5: To find the transfer function of two port network and plot the magnitude and phase plot

# Forming the transfer function

ω = 1.5 # Default Value
R = 100 
L = 1e-6 # 1e-6 <==> 10^(-6)
C = 1e-6
ω = 1/(sqrt(L * C)) # Updating omega(ω)
Q = (1/R) * (sqrt(L/C)) # Quality factor 
ζ = 1/(2 * Q) # Damping factor

num = poly1d([ω**2])
den = poly1d([1, 2*ω*ζ, ω**2])

H5 = sp.lti(num, den)

ω, mag, phi = H5.bode()

# Plotting the magnitude and phase plot

figure(4) # Magnitude Plot
semilogx(ω, mag)
title("Q5: Magnitude Bode plot")
xlabel("ω  → ")
ylabel("20log|H(jω)| → ")
grid(True)

figure(5) # Phase plot
semilogx(ω, phi)
title("Q5: Phase Bode plot")
xlabel("ω  → ")
ylabel("∠H(jω) → ")
grid(True)
show()

# Q6: To Find the output voltage from transfer function and input voltage for short term and long term time intervals, and plot for the same.

# For short time interval

t6 = arange(0, 30e-6, 1e-8)
vi = cos(1e3 * t6) - cos(1e6 * t6)
t6, vo_short, svec = sp.lsim(H5, vi, t6)

figure(6) 
plot(t6, vo_short)
title("Q6: The Output Voltage for short time interval")
xlabel("t → ", fontsize = 13)
ylabel("v\u2092(t) → ", fontsize = 13)
grid(True)

# For Long time interval

t7 = arange(0, 10e-3, 1e-8)
vi = cos(1e3 * t7) - cos(1e6 * t7)
t7, vo_long, svec = sp.lsim(H, vi, t7)

figure(7)
plot(t7, vo_long)
title("Q6: The Output Voltage for long time interval")
xlabel("t → ", fontsize = 13)
ylabel("v\u2092(t) → ", fontsize = 13)
grid(True)
show()









