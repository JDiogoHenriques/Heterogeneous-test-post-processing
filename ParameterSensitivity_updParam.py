"""
    description
    Varies a given paramter and runs the ABAQUS analysis in order to perform a sensitivity 
    analysis.

    This script is only prepared for the Hill48+Swift law constitutive model and for the
    Load90_MatDir90 configuration.
"""

# IMPORT PACKAGES -------------------------------------------------------------
from pickle import TRUE
from PythonScript_main import *
import numpy as np
import os
import time
from subprocess import Popen,PIPE,TimeoutExpired

#  --------------------------------------------------------------- method - MAIN
def main():
    ref_param_input = np.zeros(len(Ref_param)+1)
    ref_param_input[0] = Ref_param[0]
    ref_param_input[1] = 1-Ref_param[1]
    for i in range(2,len(Ref_param)+1):
        ref_param_input[i] = Ref_param[i-1]

    cwd = os.getcwd()
    cwd_sens = os.path.join(cwd,test_config,'ParameterSensitivity\\FEA\\')


    j=0
    for j in range(len(Ref_param)):
        cw_senseParam = os.path.join(cwd_sens,Param[j])
        os.chdir(cw_senseParam)
        ref_param_variat = np.copy(ref_param_input)

        if j==1:
            ref_param_variat[j+1]= ref_param_input[j+1]-(ref_param_input[j+1]*Param_variation)
            ref_param_variat[j]= ref_param_input[j]+(ref_param_input[j+1]*Param_variation)
        elif j>1:
            ref_param_variat[j+1]= ref_param_input[j+1]-(ref_param_input[j+1]*Param_variation)
        else:
            ref_param_variat[j]= ref_param_input[j]-(ref_param_input[j]*Param_variation)

        try:
            with open('Param.inp','w') as f:
                f.write('*Parameter')
                for k in range(0,len(ref_param_variat)):
                    f.write('\n')
                    f.write('V{:d}={:.12f}'.format(k+1,ref_param_variat[k]))
            print('Parameters of analysis ' + str(Param[j]) + ' was updated.')
        except:
            print('Failed to update parameters on analysis' + str(Param[j]))
            pass

    print('Done.')        

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
