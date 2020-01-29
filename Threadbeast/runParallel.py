import sys
import time
from math import *
import numpy as np
import os
import concurrent.futures
for j in range(8,16):
    valTemp = 0.7 - j*0.07
    Dn = exp(valTemp)
    Dc = Dn*0.3
    for k in range(0,4):
        Chi = 1.0 + k*0.6
        logpath = "logs/Dn" + str(j) + "Chi" + str(k)
        # print(logpath)
        command="sim/build/processNoDisplay " + logpath +" 0.00001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 500000" + " 100000 " + " 0 >output &"
        # print(command)
        #  os.system(command)
        with concurrent.futures.ProcessPoolExecutor(max_workers=32) as executor:
          executor.submit(os.system(command))
        #  executor.submit(command)
print("sript finished")

