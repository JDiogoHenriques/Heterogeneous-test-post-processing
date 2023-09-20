"""
    description
   Extraction of the Major and minor strain components

    INPUT: .odb file of the test with the LE, SDV1, S and EVOL field output request 
    OUTPUT: .csv files with principal stresses and strains, eq. plastic strain, and element volume


    Run by writting on the command line: "abaqus cae noGUI=AbaqusPythonScript_Results.py"
"""

# IMPORT PACKAGES -------------------------------------------------------------
from PythonScript_main import *
import numpy as np
from pyDOE import *
import os

#  --------------------------------------------------------------- method - MAIN
def main():
    cwd = os.getcwd()
    path_to_save_Start = (cwd + "\\" + test_config + "\\MatchID\\VFM\\StartingParam.csv")
    path_to_save = (cwd + "\\" + test_config + "\\MatchID\\VFM\\LB_UB.csv")
    number_runs = 4
    Param = Ref_param
    UB = np.round(np.array(Ref_param)*1.3, 4)
    LB = np.round(np.array(Ref_param)*0.7, 4)
    Design = lhs(len(Param), samples=number_runs, criterion='center')

    StartingParam = np.zeros((number_runs,len(Param)))
    for i in range(number_runs):
        StartingParam[i,:] = np.round((np.array(UB-LB))*Design[i,:] + np.array(LB),4)

    with open(path_to_save_Start, 'w') as f:
        np.savetxt(f, StartingParam, fmt='%.6f')

    data = np.vstack((LB, UB))

    with open(path_to_save, 'w') as f:
        np.savetxt(f, data, fmt='%.6f')
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
