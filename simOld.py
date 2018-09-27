import sys
import time
import processes as P
from math import *
import numpy as np

basePort = 8000
'''
worlds = [P.pTemp('processes/sim/build','world00','logs/log00',1,basePort+0),
          P.pTemp('processes/sim/build','world01','logs/log01',2,basePort+1]
'''
worlds = [P.pTemp('processes/sim/build','world01','logs/log01',1,basePort+1)]


sizen1 = 20.0
'''
Dn = np.array([sizen1,sizen2])
Dc = np.array([sizen1*0.3,sizen2*0.3])
'''
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
for t in range(400000):


    for i,w in enumerate(worlds):
        w.stream('1,')

    for i,w in enumerate(worlds):
        if(t%200  == 0):
            w.stream('2,') # display image

        if (t%40000==0):
            w.stream('3,') #Save image to file
            print(t)

 # Uncomment to save
for i,w in enumerate(worlds):
    w.stream('5,'+'logs/data'+str(i)+'.bin')
w.quit()
