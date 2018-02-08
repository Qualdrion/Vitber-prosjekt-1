from math import *
from matplotlib import pyplot as plt

#Function that takes in the d-value and creates "points" y-values to use for plotting using the formula for the y-values
#from the assignment file.
def Generate_Y (d, points):
    y = []
    for i in range (points+1):
        x = i/points
        y.append((2/3-x)/(d*(d**2+(x-2/3)**2)**(1/2))-(1/3-x)/(d*(d**2+(x-1/3)**2)**(1/2)))
    return y

#Function used to create the points on the x-axis. Pretty unnecessary considering there are similar methods in numpy and such
#but now that I made it for some reason there isn't really any reason to remove it.
def Generate_X (points):
    x = []
    for i in range (points+1):
        x.append(i/points)
    return x

#For x it loops through a list y and plots the points x, yi, for all yi in the list y. It uses a logarithmic scale on the
#y-axis as specified in the assignment.
def PlotXY (x, y):
    plt.plot(x, y[0], label = 'd=0.025')
    plt.plot(x, y[1], label = 'd=0.25')
    plt.plot(x, y[2], label = 'd=2.5')
    plt.legend()
    plt.semilogy()
    plt.show()

#Code to use the methods to generate the desired plot from problem 1.
d = [0.025, 0.25, 2.5]
y = []
for i in range(len(d)):
    y.append(Generate_Y(d[i], 100))
x = Generate_X(100)

PlotXY(x, y)