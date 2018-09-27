import sys
import time
import processes as P
from math import *
import numpy as np

basePort = 8000
'''
worlds = [P.pTemp('processes/simWhole/build','world00','logs/log00',1,basePort+0),
          P.pTemp('processes/simWhole/build','world01','logs/log01',2,basePort+1]
'''
worlds = [P.pTemp('processes/simWhole/build','world01','logsWhole/log01',1,basePort+1)]

sizen1 = 60
'''
Dn = np.array([sizen1,sizen2])
Dc = np.array([sizen1*0.3,sizen2*0.3])
'''
Dn = np.array([8])
Dc = np.array([2])

# Set simulation params
for i,w in enumerate(worlds):
    w.stream('4,0,'+str(Dn[i]))
    w.stream('4,1,'+str(Dc[i]))

# Uncomment to load
'''
for i,w in enumerate(worlds):
    w.stream('6,'+'logsWhole/data'+str(i)+'.bin')
'''
# MAIN SIMULATION LOOP

for t in range(2000):

    for i,w in enumerate(worlds):
        w.stream('1,')

    for i,w in enumerate(worlds):
        if(t%20  == 0):
            w.stream('2,') # display image
        if (t%1000==0):
            w.stream('3,') #Save image to file
            print(t)

 # Uncomment to save
for i,w in enumerate(worlds):
    w.stream('5,'+'logsWhole/data'+str(i)+'.bin')
w.quit()
