#!/usr/bin/python3

import os

pathname = os.path.abspath('.')
savePath = os.path.join(pathname,'postLine')
if not os.path.isdir(savePath):
  os.makedirs(savePath)

postPath = os.path.join(pathname,'postProcessing/sets')
if not os.path.isdir(postPath):
  postPath = 'postProcessing/line'

# List of time dirs in order
a = sorted(os.listdir(postPath), key=float)

# Get number of sensors
dir1 = os.path.join(pathname,postPath,a[int(len(a)/2.0)]) # Choose one directory
b = sorted(os.listdir(dir1), key=float) # In this directory, we have sensor data
nSens = 0
index = []
for i in range(len(b)):
  test1 = b[i].find('line') + 1
  test2 = b[i].find('U') + 1
  if test1 and test2:
    index.append(i)
    nSens += 1

first = True

for i in range(nSens):
  # Create files to write
  #fileName = b[index[i]][0:b[index[i]].find('_')] # e.g., line1
  fileName = 'probeU' + str(i) 
  fileW = open(os.path.join(savePath,fileName), 'w')
  print('Sensor ' + '%i' % int(i+1) + ' of ' + '%i' % nSens + '.')

  # Read files time by time
  for j in range(len(a)):
    directory = os.path.join(pathname,postPath,a[j])
    try:
      fileR = open(os.path.join(directory,b[index[i]]), 'r')
    except:
      print('WARNING - File not present: ' + os.path.join(directory,b[index[i]]))
    else:
      data = fileR.read()
      fileR.close()
      data = data.split('\n')

      if first: # First time step
        coord = j
        first = False

#      y = []
#      ux = []
#      uy = []
#      uz = []

      # x y z U caculation
      for k in range(len(data)-1):
        line = data[k]
        line = line.split('\t')
        
#        y.append(float(line[0]))
#        ux.append(float(line[1]))
#        uy.append(float(line[2]))
#        uz.append(float(line[3]))
        
        tmpy  = float(line[0])
        tmpux = float(line[1])
        tmpuy = float(line[2])
        tmpuz = float(line[3])

        time = a[j]
        fileW.write(time + '\t' + '%.6f' %tmpy + '\t' + '%.6f' %tmpux + '\t' + '%.6f' %tmpuy + '\t' + '%.6f' %tmpuz + '\n')

#      # Integrate in Z
#      wLevel = y[0] # Starting locatino of the probe
#      for k in range(len(y)-1):
#        wLevel = wLevel + alpha[k]*(y[k+1]-y[k]) # Write to file
#      time = a[j]
#      fileW.write(time + ' ' + '%.6f' % wLevel + '\n')

  fileW.close()

print('Done')
