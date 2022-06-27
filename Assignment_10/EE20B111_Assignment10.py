"""
            EE2703 Applied Programming Lab - 2022
            Assignment 10: Linear And Circular Convolution
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 27-04-2022
"""

# Importing the required modules
import numpy as np
import scipy.signal as sp
from pylab import *
import csv


# Q1: Reading the contents in csv file  
# h.csv contains the FIR Filter coeffs, and x1.csv contains the zadoff-chu sequence
file1 = "h.csv"

b = np.zeros(12)
i = 0
with open(file1, 'r') as f1:  
    for line in f1:
        b[i] = float(line) # Typecasting it into float type of variable
        i += 1 # To hold the count of number of values in file


# Function that plots the graph but doesn't display it( i.e., show is not used)
def make_graph(xVal, yVal, xL = 'ω', yL = 'Magnitude', head = 'xyz', save = 'x.png'):
    plot(xVal, yVal)
    xlabel(xL)
    ylabel(yL)
    title(head)
    grid(True)
    savefig(save) 
    

# Q-2: Magnitude and phase response of low pass filter 
w, h = sp.freqz(b) # Here, a = 1
# Computing the frequency response of the FIR filter at various "w" (nnormalized range of "w" = [0, π)
# H(e^jw) = B(e^jw)/A(e^jw) = b[0] + b[1]e^(-jw) + b[2]e^(-2jw) + ... + b[M]e^(-jwM) # A = 1
figure(0)
subplot(2, 1, 1)
make_graph(w, abs(h), '', '|H| (Magnitude) →', 'Q2: Magnitude and phase response for Low pass filter')
subplot(2, 1, 2)
make_graph(w, angle(h), 'ω →', '∠H (phase) →' , '')
savefig("Ass10_Figure_1.png"), show()


# Q-3: Plotting the given sequence x(.)
n = array(range(2**10))
x = cos(0.2*pi*n) + cos(0.85*pi*n) # Generating the signal 
make_graph(n, x, 'n →', 'x →', 'Q3: Plot of sequence, x = cos(0.2\u03C0n) + cos(0.85\u03C0n)', "Ass10_Figure_2.png"), show()


# Q-4: Plotting y(.) := x(.)*b(.) {convolution of x and b}
y = np.zeros(len(x))
# Loop for convolution
for i in arange(len(x)):
    for k in arange(len(b)):
        y[i] += x[i-k]*b[k]
    
make_graph(n, y, 'n →', 'y →', 'Q4: Output of linear convolution: y(n) = x(n) * b(n)', "Ass10_Figure_3.png"), show()


# Q-5: Output of circular convolution
y = ifft(fft(x)*fft(concatenate((b, zeros(len(x) - len(b))))))
make_graph(n, real(y), 'n →', 'Re{y} →', 'Q5: Output of circular convolution', "Ass10_Figure_4.png"), show()


# Q-6: Output of circular convolution using liner convolution
def circular_conv(x, h):
    P = len(h)
    n_temp = int(ceil(log2(P)))
    h_temp = np.concatenate((h, np.zeros(int(2**n_temp) - P)))
    P = len(h_temp)
    n1 = int(ceil(len(x)/2**n_temp))
    x_temp = np.concatenate((x, np.zeros(n1*(int(2**n_temp)) - len(x))))
    y = np.zeros(len(x_temp) + len(h_temp) - 1)
    for i in range(n1):
        temp = np.concatenate((x_temp[i*P:(i + 1)*P], np.zeros(P - 1)))
        y[i*P:(i + 1)*P + P - 1] += np.fft.ifft(np.fft.fft(temp) * np.fft.fft( np.concatenate( (h_temp,np.zeros(len(temp)-len(h_temp))) ))).real
    return y

y = circular_conv(x, b)
make_graph(n, real(y[:1024]), 'n →', 'Re{y} →', 'Q6: Output of circular convolution using linear convolution', "Ass10_Figure_5.png"), show()


# Q-7: AUTO-CORRELATION
file2 = "x1.csv"

lines1 = []
with open(file2, 'r') as f2:
    csvreader = csv.reader(f2)
    for row in csvreader:
        lines1.append(row)

lines2 = []
for line in lines1:
    line = list(line[0])
    try :
        line[line.index('i')] = 'j'
        lines2.append(line)
    except ValueError:
        lines2.append(line)
        continue
x = [complex(''.join(line)) for line in lines2]
X = np.fft.fft(x)
x2 = np.roll(x, 5)
cor = np.fft.ifftshift(np.correlate(x2, x, 'full'))
print("The length of correlation array of x1 and shifted version of x1: ", len(cor))

figure()
xlim(0, 20)
make_graph(linspace(0, len(cor) - 1, len(cor)), abs(cor), 't →', 'Correlation →', 'Q7: Auto-Correlation of x1 and shifted version(right shift by 5) of x1', "Ass10_Figure_6.png"), show()
