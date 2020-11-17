import sys
import time
from math import *
import numpy as np
import os
for j in range(0,99):
    valTemp = 3.6
    Dn =  exp(valTemp)
    Dc = 0.3*Dn
    for k in range(1,2):
        Chi = (0.0 + k*0.4)*Dn
        seed = j  * 1389 + 1
        logpath = "./logsStaticDn/Dn" + str(j);
        os.system("mkdir " + logpath);
        for m in range(0,1):
            reg = str(m)
            #os.system("cp logsStaticDn/13.h5 " + logpath)
            '''
            Lcontinue LfixedSeed Lgraphics LDn
            '''
            command="$HOME/Neuroscience/DirichletRD/Threadbeast/sim/build/processSingle " + logpath + " " + reg + " 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 8 " + " 5.0 " +  " 200000" + " 199999 " + " 0 " + str(seed) + " 1 " + " 0  > " + logpath + "/output "
            print(logpath)
            print(command)
            os.system(command)
            time.sleep(10)
print("script finished");
