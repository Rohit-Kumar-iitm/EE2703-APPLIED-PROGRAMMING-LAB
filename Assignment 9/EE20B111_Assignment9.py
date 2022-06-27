"""
            EE2703 Applied Programming Lab - 2022
            Assignment 9: Spectra Of Non-Periodic Signals
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 16-04-2022
"""

# Importing the required modulues
from pylab import*
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# Assignment Question 1

def Example1():
    t = linspace(-pi, pi, 65)[:-1]
    dt = t[1] - t[0] 
    fmax = 1/dt
    y = sin(sqrt(2)*t)
    y[0] = 0 # The sample corresponding to -tmax should be set to zero
    y = fftshift(y) # Make y start with y(t=0)
    Y = fftshift(fft(y))/64.0
    w = linspace(-pi*fmax, pi*fmax, 65)[:-1]
    
    figure()
    subplot(2, 1, 1)
    title("Spectrum of sin(\u221A2t)") # utf code for square root -> U+221A
    ylabel("|Y|", size = 15)
    xlim([-10, 10])
    plot(w, abs(Y), lw = 2)
    grid(True)
    subplot(2, 1, 2)
    xlabel("\u03C9", size = 15)
    ylabel("Phase of Y", size = 15)
    xlim([-10, 10])
    plot(w, angle(Y), 'ro', lw = 2)
    grid(True)
    savefig("fig10-1.png")
    show() 
    
def Example2():
    t1 = linspace(-pi, pi, 65)[:-1]
    t2 = linspace(-3*pi, -pi, 65)[:-1]
    t3 = linspace(pi, 3*pi, 65)[:-1]
    # y = sin(sqrt(2)*t)
    
    figure(2)
    title("sin(\u221A2t)")
    xlabel("t", size = 15)
    ylabel("y", size = 15)
    plot(t1, sin(sqrt(2)*t1), 'b', lw = 2)
    plot(t2, sin(sqrt(2)*t2), 'r', lw = 2)
    plot(t3, sin(sqrt(2)*t3), 'r', lw = 2)
    grid(True)
    savefig("fig10-2.png")
    show()

def Example3():
    t1 = linspace(-pi, pi, 65)[:-1]
    t2 = linspace(-3*pi, -pi, 65)[:-1]
    t3 = linspace(pi, 3*pi, 65)[:-1]
    y = sin(sqrt(2)*t1)
    
    figure(3)
    title("sin(\u221A2t) with t wrapping every 2\u03C0")
    xlabel("t", size = 15)
    ylabel("y", size = 15)
    plot(t1, y, 'bo', lw = 2)
    plot(t2, y, 'ro', lw = 2)
    plot(t3, y, 'ro', lw = 2)
    grid(True)
    savefig("fig10-3.png")
    show()

def GibbsPhenomenon():
    t = linspace(-pi, pi, 65)[:-1]
    dt = t[1] - t[0]
    fmax = 1/dt
    y = t
    y[0] = 0 # The sample corresponding to -tmax should be set to zero
    y = fftshift(y) # Make y start with y(t=0)
    Y = fftshift(fft(y))/64.0
    w = linspace(-pi*fmax, pi*fmax, 65)[:-1]
    
    figure()
    title("Spectrum of a digital ramp")
    xlabel("\u03C9", size = 15)
    ylabel("|Y| (in dB)", size = 15)
    xlim([1, 10])
    ylim([-20, 0])
    semilogx(abs(w), 20*log10(abs(Y)), lw = 2)
    xticks([1, 2, 5, 10],["1", "2", "5", "10"], size = 16)
    grid(True)
    savefig("fig10-4.png")
    show() 

def Windowing1():
    t1 = linspace(-pi, pi, 65)[:-1]
    t2 = linspace(-3*pi, -pi, 65)[:-1]
    t3 = linspace(pi, 3*pi, 65)[:-1]
    k = arange(64)
    wnd = fftshift(0.54 + 0.46*cos(2*pi*k/63))
    y = sin(sqrt(2)*t1)*wnd
    
    figure(3)
    title("sin(\u221A2t)xw(t) with t wrapping every 2\u03C0")
    xlabel("t", size = 15)
    ylabel("y", size = 15)
    plot(t1, y, 'bo', lw = 2)
    plot(t2, y, 'ro', lw = 2)
    plot(t3, y, 'ro', lw = 2)
    grid(True)
    savefig("fig10-5.png")
    show()   


def Windowing2():
    t = linspace(-pi, pi, 65)[:-1]
    dt = t[1] - t[0]
    fmax = 1/dt
    k = arange(64)
    wnd = fftshift(0.54 + 0.46*cos(2*pi*k/63))
    y = sin(sqrt(2)*t)*wnd
    y[0] = 0 # The sample corresponding to -tmax should be set to zero
    y = fftshift(y) # Make y start with y(t=0)
    Y = fftshift(fft(y))/64.0
    w = linspace(-pi*fmax, pi*fmax, 65)[:-1]
    
    figure()
    subplot(2, 1, 1)
    title("Spectrum of sin(\u221A2t) x w(t)")
    ylabel("|Y|", size = 15)
    xlim([-8, 8])
    plot(w, abs(Y), lw = 2)
    grid(True)
    subplot(2, 1, 2)
    xlabel("\u03C9",size = 15)
    ylabel("Phase of Y", size = 15)
    xlim([-8, 8])
    plot(w, angle(Y), 'ro', lw = 2)
    grid(True)
    savefig("fig10-6.png")
    show()     

def Windowing3():
    t = linspace(-4*pi, 4*pi, 257)[:-1]
    dt = t[1] - t[0]
    fmax = 1/dt
    k = arange(256)
    wnd = fftshift(0.54 + 0.46*cos(2*pi*k/256))
    y = sin(sqrt(2)*t)
    y = y*wnd
    y[0] = 0 # the sample corresponding to -tmax should be set zeroo
    y = fftshift(y) # make y start with y(t=0)
    Y = fftshift(fft(y))/256.0
    w = linspace(-pi*fmax, pi*fmax, 257)[:-1]
    figure()
    subplot(2, 1, 1)
    title("Spectrum of sin(\u221A2t) x w(t)")
    ylabel("|Y|$", size = 15)
    xlim([-8, 8])
    plot(w, abs(Y), lw = 2)
    grid(True)
    subplot(2, 1, 2)
    xlabel("\u03C9",size = 15)
    ylabel("Phase of Y", size = 15)
    xlim([-8, 8])
    plot(w, angle(Y), 'ro', lw = 2)
    grid(True)
    savefig("fig10-7.png")
    show()          

# Plotting graphs
Example1()
Example2()
Example3()
GibbsPhenomenon()
Windowing1()
Windowing2()
Windowing3()

# Helper function
def spectrum(lim, n, f, t_temp = 0, show_temp = True, t_lims = False, windowing = False, xlim1 = 10, title1 = "Spectrum of sin(\u221A2t)", xlabel1 = "\u03C9", ylabel1 = "|Y|", ylabel2 = "Phase of Y", savename = "abc.png"):
    if(t_lims):
        t = t_temp
    else:
        t = linspace(-lim, lim, n + 1)[:-1]
    dt = t[1] - t[0]
    fmax = 1/dt
    y = f(t)
    if(windowing):
        m = arange(n)
        wnd = fftshift(0.54 + 0.46*cos(2*pi*m/n))
        y = y*wnd
    y[0] = 0 # The sample corresponding to -tmax should be set to zero
    y = fftshift(y) # Make y start with y(t=0)
    Y = fftshift(fft(y))/float(n)
    w = linspace(-pi*fmax, pi*fmax, n + 1)[:-1]
    
    mag = abs(Y)
    phase = angle(Y)
    if(show_temp):
        figure()
        subplot(2, 1, 1)
        title(title1)
        ylabel(ylabel1, size = 15)
        xlim([-xlim1, xlim1])
        plot(w, mag, lw = 2)
        grid(True)
        subplot(2, 1, 2)
        phase[where(mag < 3e-3)] = 0
        xlabel(xlabel1, size = 15)
        ylabel(ylabel2, size = 15)
        xlim([-xlim1, xlim1])
        plot(w, phase, 'ro', lw = 2)
        grid(True)
        savefig(savename)
        show()
    return w, Y

# Defining cos^3(t)
def cosCube(t, w0 = 0.86):
    return (cos(w0*t))**3

# Assignment Question 2

# FFT of cos^3(t) windowed and unwindowed
a, b = spectrum(4*pi, 64*4, cosCube, xlim1 = 3, windowing = False, title1 = "Spectrum of cos\u00B3(\u03C9\u2092t)", savename = '10-8_NoWindow.png')
a, b = spectrum(4*pi, 64*4, cosCube, xlim1 = 3, windowing = True, title1 = "Spectrum of cos\u00B3(\u03C9\u2092t)", savename = '10-8_Window.png')

def cosine(t, w0 = 1.5, delta = 0.5):
    return cos(w0*t + delta)

# Assignment Question 3
#FFT of cos(wt + delta) windowed to estimate w, delta
w, Y = spectrum(pi, 128, cosine, xlim1 = 3, windowing = True, title1 = "Spectrum of cos(\u03C9\u2092t + \u03B4)", savename = '10-8.png')

def estimated_omega(w, Y):
    ii = where(w > 0)
    omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2)) # Weighted average
    print ("Omega = " + str(omega))

def estimated_delta(w, Y, sup = 1e-4, window = 1):
    ii_1 = np.where(np.logical_and(np.abs(Y) > sup, w > 0))[0]
    np.sort(ii_1)
    points = ii_1[1:window + 1]
    print("Delta = " + str(np.sum(np.angle(Y[points]))/len(points))) # Weighted average for first 2 points

print("The estimated value of \u03C9 and \u03B4 from FFT of cos(\u03C9t + \u03B4) are: ")
estimated_omega(w, Y)
estimated_delta(w, Y)

# Assignment Question 4
def noisecosine(t, w0 = 1.5, delta = 0.5):
    return cos(w0*t + delta) + 0.1*np.random.randn(128) # Adding noise (treating noise as a random variable)

w, Y = spectrum(pi, 128, noisecosine, xlim1 = 3, windowing = True, title1 = "Spectrum of cos(\u03C9\u2092t + \u03B4)")

print("The estimated value of \u03C9 and \u03B4 from FFT of cos(\u03C9t + \u03B4) + noise are: ")
estimated_omega(w, Y)
estimated_delta(w, Y)

# Assignment Question 5
 
def chirp(t):
    return cos(16*(1.5 + t/(2*pi))*t) 

w, Y = spectrum(pi, 1024, chirp, xlim1 = 60, windowing = True, title1 = "Spectrum of chirp function")
w, Y = spectrum(pi, 1024, chirp, xlim1 = 60, windowing = False, title1 = "Spectrum of chirp function")

# Assignmnt Question 6

t = np.linspace(-np.pi, np.pi, 1025); t = t[:-1]
t_arrays = np.split(t, 16)

Y_mags = np.zeros((16, 64))
Y_angles = np.zeros((16, 64))

# Splitting array and performing FFT
for i in range(len(t_arrays)):
    w, Y = spectrum(lim = 10, t_temp= t_arrays[i], t_lims = True, n = 64, f = chirp, xlim1 = 60, windowing = False, title1 = "Spectrum of chirp function", show_temp = False)
    Y_mags[i] = abs(Y)
    Y_angles[i] = angle(Y)
# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

t = np.linspace(-np.pi, np.pi, 1025); t = t[:-1]
fmax = 1/(t[1] - t[0])
t = t[::64]
w = np.linspace(-fmax*np.pi, fmax*np.pi, 64 + 1); w = w[:-1]
t,w = np.meshgrid(t, w)

surf = ax.plot_surface(w, t, Y_mags.T, cmap = cm.coolwarm, linewidth = 0, antialiased = False)
fig.colorbar(surf, shrink = 0.5, aspect = 5)
plt.ylabel("Frequency")
plt.xlabel("Time")

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
surf = ax.plot_surface(w, t, Y_angles.T, cmap = cm.coolwarm, linewidth = 0, antialiased = False)
fig.colorbar(surf, shrink = 0.5, aspect = 5)
plt.ylabel("Frequency")
plt.xlabel("Time")

plt.show()