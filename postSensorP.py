#!/usr/bin/python3

import os

pathname = os.path.abspath('.')
savePath = os.path.join(pathname,'postLine')
if not os.path.isdir(savePath):
  os.makedirs(savePath)

postPath = os.path.join(pathname,'postProcessing/sets')
if not os.path.isdir(postPath):
  postPath = 'postProcessing/GaugesP' # Directory name for pressure

# List of time dirs in order
a = sorted(os.listdir(postPath), key=float)

# Get number of sensors
dir1 = os.path.join(pathname,postPath,a[int(len(a)/2.0)]) # Choose one directory
b = sorted(os.listdir(dir1), key=float) # In this directory, we have sensor data
nSens = 0
index = []
for i in range(len(b)):
  test1 = b[i].find('Gauges') + 1
  test2 = b[i].find('p') + 1
  if test1 and test2:
    index.append(i)
    nSens += 1

first = True

for i in range(nSens):
  # Create files to write
  #fileName = b[index[i]][0:b[index[i]].find('_')] # e.g., line1
  fileName = 'gauge' + str(i+1) 
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
      nProbe = len(data)-1 # Number of gagues

      if first: # First time step
        # Save the probe location
        for p in range(0,3):
          for k in range(nProbe):
            if ( k == 0 ):
              fileW.write('0' + '\t') # dummy zeros for time column
            pLoc = data[k]
            pLoc = pLoc.split('\t')
            fileW.write(pLoc[p]+'\t')
            if (( (k+1) % nProbe ) == 0 ):
              fileW.write('\n')
        
       # coord = j
        first = False


      # x y z alpha1 caculation
      for k in range(len(data)-1):
        if ( k == 0 ):
          time = a[j]
          fileW.write(time + '\t')
        line = data[k]
        line = line.split('\t')
        fileW.write(line[3]+'\t')
        if (( (k+1) % nProbe ) == 0 ):
          fileW.write('\n')

  fileW.close()

print('Done')
