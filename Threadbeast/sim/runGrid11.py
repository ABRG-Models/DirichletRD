import sys
import time
from math import *
import numpy as np
import os
for j in range(0,8):
    valTemp = 3.6 -  j*0.4
    Dn = exp(valTemp)
    Dc = 0.3*Dn
    for k in range(1,2):
        Chi = (0.0 + k*0.4)*Dn
        logpath = "./logs11" + str(j);
        os.system("mkdir " + logpath);
        os.system("cp ./logs/11.h5 " + logpath);
        command="$HOME/Neuroscience/DirichletRD/Threadbeast/sim/build/processSingle " + logpath + " 11 " + " 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 8 " + " 5.0 " +  " 400000" + " 399999 " + " 0 " + " 0 " + " 1 " + " 0  > " + logpath + "/output "
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
