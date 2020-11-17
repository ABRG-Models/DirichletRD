import sys
import time
from math import *
import numpy as np
import os
for j in range(0,8):
    valTemp = 3.6 -  j*0.0125
    Dn = exp(valTemp)
    Dc = 0.3*Dn
    for k in range(1,2):
        Chi = (0.0 + k*0.4)*Dn
        logpath = "./logsStaticDn/Dn" + str(j);
        os.system("mkdir " + logpath);
        command="$HOME/Neuroscience/DirichletRD/Threadbeast/sim/build/processSetup " + " ./logsStaticDn " +  " 8 " + " 5.0 " + " 0  1 > output "
        print(logpath)
        print(command)
        os.system(command)
        for m in range(0,21):
            os.system("cp ./logsFine/" + str(m) + ".h5 " + logpath);
print("script finished");
