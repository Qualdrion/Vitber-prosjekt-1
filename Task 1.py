from math import *
from matplotlib import pyplot as plt

def Generate_Y (d, points):
    y = []
    for i in range (points+1):
        x = i/points
        y.append((2/3-x)/(d*(d**2+(x-2/3)**2)**(1/2))-(1/3-x)/(d*(d**2+(x-1/3)**2)**(1/2)))
    return y

def Generate_X (points):
    x = []
    for i in range (points+1):
        x.append(i/points)
    return x

def PlotXY (x, y):
    for i in range (len(y)):
        plt.plot(x, y[i])
    plt.semilogy()
    plt.show()

d = [0.025, 0.25, 2.5]
y = []
for i in range(len(d)):
    y.append(Generate_Y(d[i], 100))
x = Generate_X(100)

PlotXY(x, y)