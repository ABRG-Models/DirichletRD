import sys
import time
from math import *
import numpy as np
import os
for j in range(1,3):
    valTemp = 3.6 -  j*0.05
    Dn = exp(valTemp)
    Dc = 0.3*Dn
    for k in range(0,1):
        Chi = (0.0 + k*0.4)*Dn
        logpath = "./logs9" + str(j);
        os.system("mkdir " + logpath);
        os.system("cp ./logs/9.h5 " + logpath);
        command="$HOME/Neuroscience/DirichletRD/Threadbeast/sim/build/processField " + logpath + "9" + " 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 8 " + " 5.0 " +  " 100000" + " 99999 " + " 0 " + " 1 " + " 1 " + " 0  > " + logpath + "/output "
        print(logpath)
        print(command)
        os.system(command)
        time.sleep(10)
        '''
    for k in range(10,32):
        Chi = (1.0 + k*0.064)*Dn
        logpath = "./Ratlogs/Random064/Dn" + str(j) + "Chi" + str(k)
        command="./processField " + logpath +" 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 8 " + " 5.0 " +  " 100" + " 95 " + " 1 > " + logpath + "/output "
        print(logpath)
        print(command)
        os.system(command)
        time.sleep(10)
        '''
print("script finished");
