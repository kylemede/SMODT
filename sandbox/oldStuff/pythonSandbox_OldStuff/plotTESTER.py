import math
import numpy as np
import os
import pylab
from numpy import linspace
from scipy import pi

plt = pylab.matplotlib.pyplot
mpl = pylab.matplotlib.mpl

import random as rand

zLabel='zlabel'
xLim=False
yLim=False
xLabel='xLabel'
yLabel='yLabel'

xData=[]
yData=[]
zData=[]
xData.append(1)
yData.append(1)
zData.append(0.001)
for i in range(0,1000):
    xData.append(xData[-1]+0.01)
    yData.append(rand.uniform(1,200))
    #zData.append(zData[-1]+0.001)
    zData.append(zData[-1]+0.1)





