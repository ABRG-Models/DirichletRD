import sys
import time
import processes as P
from math import *
import numpy as np
for j in range(0,10):
  basePort = 8000
  worlds = [P.pTemp('processes/sim/build',str(j),'logs/'+str(j),1,basePort+1)]
  valTemp = 3.0 - j*0.1

  sizen1 = exp(valTemp)
  Dn = np.array([sizen1])
  Dc = np.array([0.3*sizen1])

# Set simulation params
  for i,w in enumerate(worlds):
      w.stream('4,0,'+str(Dn[i]))
      w.stream('4,1,'+str(Dc[i]))

# Uncomment to load
  for i,w in enumerate(worlds):
      w.stream('6,'+'logs/data'+str(i)+'.bin')

# MAIN SIMULATION LOOP
  for t in range(50000):


      for i,w in enumerate(worlds):
          w.stream('1,')

      for i,w in enumerate(worlds):
          if(t%200  == 0):
              w.stream('2,') # display image

          if (t%49998==0):
              w.stream('3,') #Save image to file
              print(t)

 # Uncomment to save
  for i,w in enumerate(worlds):
      w.stream('5,'+'logs/data'+str(i)+'.bin')
  w.quit()
