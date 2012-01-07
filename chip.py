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
