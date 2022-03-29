#!/usr/bin/python3

import os 
from pylab import *

pathname = os.path.abspath('.')
readPath = os.path.join(pathname,'postLine')
a = os.listdir(readPath)

# Sorting
remove = []
for i in range(len(a)):
  if (a[i].rfind('.')+1): # Includes point
     remove.append(i)
remove.reverse() # to remove the elements from behind
for i in remove:
  a.pop(i)
a = sorted(a, key=lambda v: int(v.split('e')[1])) # sort by probe numbers

# Plot  
index = 0
indexFig = 0
for gauge in a: # a is now sorted list
  index = index + 1
  if index >= 5 or index == 1: # Now, working with 4 probes
    index = 1
    indexFig = indexFig + 1
    figure(num=indexFig)

  subplots_adjust(hspace=0.6)
    
  fileR = open(os.path.join(readPath,gauge),'r')
  data = fileR.read()
  fileR.close()
  data = data.split('\n')
  x = []
  y = []
  for i in range(len(data)-1):
    line = data[i]
    line = line.split(' ') 
    x.append(float(line[0]))
    y.append(float(line[1]))
              
    subplot(4,1,index)
    plot(x,y)
    xlabel('t [s]')
    ylabel('$h + \eta$ [m]')
    title(gauge)

# show()
#savefig('test.png')
