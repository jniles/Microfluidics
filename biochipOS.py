"""
A suite of tools for managing digital microfluidic biochips.
Developed by Jonathan Niles
January 2012
"""

import os
import sys
import numpy
import pulp

#####################################
# TODO: Impliment Graphical User touchpad using PyQt4
# FIXME: Find out what causes the PyQt4 error for pymodules/python2.6/sip.o
# TODO: Impliment VisChip in Chaco!
# TODO: Impliment http://code.enthought.com/projects/chaco/gallery.php selectable colormap as plot
# FIXME: Abs() line 177-178 
# 
# TODO: Integrate Droplet+Sink+Chip
#
#
#####################################

class Droplet(object):
  def __init__(self, x, y):
    """
    Creates a droplet object to navigate the
    biochip.
    """
    self.x = x
    self.y = y
    self.t = 0

  def update(self, x=None, y=None):
    """
    updates the droplet coordinates
    """
    if x:
      self.x = x
    elif y:
      self.y = y
    self.t += 1

class Sink(object):
  def __init__(self, x, y):
    """
    creates a sink on the biochip
    """
    self.x = x
    self.y = y

class Chip:
  def __init__(self, xsize, ysize, tsize):
    """
    A biochip object of dimensions x, y, and t.
    """
    self.chip = numpy.ndarray((xsize, ysize, tsize),                    #creates an array with values and
                dtype=[('v',int),('t',int)])                            #timestamps

    for dtype in self.chip.dtype.descr:
      for x in range(self.chip.shape[0]):
        for y in range(self.chip.shape[1]):
          for t in range(self.chip.shape[2]):
            if dtype[0] == 'v':
              self.chip[dtype[0]][x][y][t] = 0                          #sets all values in the array to 0
            else:
              self.chip[dtype[0]][x][y][t] = t+1                        #sets time values in the array 0=<i=<t
  
  def timeslice(self,t):
    """
    a flat slice of the array at time t
    """
    if t > self.chip.shape[2] or t < 1:
      raise NameError
    flat = numpy.zeros(self.chip.shape[:-1], dtype=int)                 #shape is identical
    for x in range(self.chip.shape[0]):
      flat[x] = [self.chip['v'][x][y][t-1] for y in range(self.chip.shape[1])] #sets flat rows equal to x rows at time t
    return flat

class Concurrency:
  def __init__(self, biochip, source, sink, path):
    """
    initializes a flat representation of the chip,
    constructs the ILP model and awaits execution.
    biochip =>  The 3D biochip instance.
    source => A tuple (a,b) of the location of the source
    sink => A tuple (i,j) of the location of the sink
    path => A list of tuples [(a,b), (c,d) ,..., (n,m)]
            that denote off-limit cells.
    """
    x,y = biochip.chip.shape[:-1]
    t = len(path)
    m,n = sink
    xst,yst = source
    #Create biochip as array of dictionaries
    chip = pulp.LpVariable.dicts("ArraySpot",
                      [(i,j,k) for i in range(1,x+1)
                      for j in range(1,y+1)
                      for k in range(1,t+1)],
                      0,1, pulp.LpBinary)
    #Add sink
    for i in range(1,t+1):
      chip[(n,m,i)] = pulp.LpVariable("Sink", pulp.LpBinary)

    #Add virtual cell
    for i in range(1, t+1):
      chip[(n+1,m, i)] = pulp.LpVariable("VirtualSpot", pulp.LpBinary)

    #Create problem instance
    self.ilp = pulp.LpProblem('Concurrency Testing ILP', pulp.LpMinimize)

    #OBJECTIVE FUNCTION
    self.ilp += pulp.lpSum(chip[m,n,k]*k for k in range(1,t+1))

    #CONSTRAINTS -- Testing Requirement
    #Any cell (i,j) in the array available for testing
    #should be visited at least once
    self.ilp += pulp.lpSum(chip[i,j,k] for i in range(1,x+1)
                                       for j in range(1,y+1)
                                       for k in range(1,t+1)) >= 1
 
    #Any cell (i,j) in the array that is running biomedical
    #assays cannot be visited by the test droplet
    self.ilp += pulp.lpSum(chip[i,j,k] for i in range(1,x+1)
                                       for j in range(1,y+1)
                                       for k in range(1,t+1)
                                       if (i,j) in path) == 0
 
    #The sink (n,m) should be visited by the test droplet
    #exactly once
    self.ilp += pulp.lpSum(chip[m,n,k] for k in range(1,t+1)) == 1

    #CONSTRAINTS -- Resource Constraint
    #
    self.ilp += pulp.lpSum(chip[i,j,k] for i in range(1, m+1)
                                   for j in range(1,n+1)
                                   for k in range(1,t+1)
                                   if (i,j,k) in chip.keys()) == 1

    #
    for k in range(2, t+1):
      self.ilp += chip[m+1,n,k] == pulp.lpSum(chip[m,n,x] for x in range(1, k-1))

    #CONSTRAINTS -- Starting Point
    #Source should be visited by the test droplet at t == 1
    self.ilp += chip[yst+1,yst,1] == 1

    #CONSTRAINTS == Movement Rules
    #
    rowSum = dict()
    for x in range(1, t+1):
      rowSum[x] = pulp.lpSum(chip[i,j,x]*i for i in range(1,m)
                                      for j in range(1,n))

    colSum = dict()
    for x in range(1, t+1):
      colSum[x] = pulp.lpSum(chip[i,j,k]*j for i in range(1,m)
                                      for j in range(1,n))

    deltaColSum = dict()
    deltaRowSum = dict()
    for x in range(1,t):
      deltaColSum[x] = colSum[x+1] - colSum[x] #FIXME: Absolute Value
      deltaRowSum[x] = rowSum[x+1] - rowSum[x] #FIXME: Absolute Value

    self.ilp += deltaColSum[x] + deltaRowSum[x] <= 1
    return

  def absolute(self, number):
    """
    An attempted workaround for the pulp/abs() clash
    """
    if number < 0:
      return -number
    else:
      return number
  
  def solve(self):
    """
    solves the ILP problem
    """
    return self.ilp.solve()


class VisChip:
  def __init__(self):
    """
    initializes the class of printer
    """
    pass

  def timeSlice(self, chipSlice, t):
    """
    visualizes the slice with all features
    at the time instance passed to the drawer.
    """
    xlist = []
    ylist = []
    ax = pylab.gca()
    for y, row in enumerate(chipSlice):
      for x, i in enumerate(row):
        if i == 1:
          xlist.append(x)
          ylist.append(chipSlice.shape[1]-y)
    ax.scatter(xlist,ylist, s=500, c='pink', marker='s')
    ax.set_xlim(0,chipSlice.shape[0])
    ax.set_ylim(0,chipSlice.shape[1])
    ax.grid()
    pylab.title('MicroFluidic BioChip Features @ T=%i'%(t))
    ax.set_xlabel('X-Coordinates')
    ax.set_ylabel('Y-Coordinates')
    pylab.show()
    return 
