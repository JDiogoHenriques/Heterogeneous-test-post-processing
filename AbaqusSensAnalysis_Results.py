"""
    Run by writting on the command line: "abaqus cae noGUI=AbaqusSensAnalysis_Results.py"
"""

# IMPORT PACKAGES -------------------------------------------------------------
from PythonScript_main import *
from abaqusConstants import *
from odbAccess import *
import numpy as np
import sys
import os
import csv


#  --------------------------------------------------------------- method - MAIN
def main():
    #Reference FEA directories
    cwd = os.getcwd()
    cwd2_ref = (test_config + '\\FEA')
    cwd_test_ref = os.path.join(cwd,cwd2_ref)

    #Perturbation FEA directories
    for n in range(len(Param)):
    #for n in range(6,6):
        cwd2 = (test_config + '\\ParameterSensitivity\\FEA\\' + Param[n])
        cwd_test = os.path.join(cwd,cwd2)
        path_to_save_force = (cwd + '\\Results\\Sensitivity_analysis\\' + test_config + '_'+Param[n]+'_Perturbation_Load.csv')
        path_and_file_results = (cwd + '\\Results\\Sensitivity_analysis\\' + test_config + '_'+Param[n]+'_Perturbation.csv')

        # ----------open odb with results of the test (it is required the COORD fieldoutput request)

        os.chdir(cwd_test)
        odb_name = test_name + '_sens.odb'
        odb = openOdb(odb_name, readOnly=True)

        errors = odb.diagnosticData.numberOfAnalysisErrors
        if errors != 0:
            return

        # ------------------------------------------------ instance of analysis
        inst = odb.rootAssembly.instances['PART-1-1']
        nset = odb.rootAssembly.instances['PART-1-1'].nodeSets['RIGHTBOUNDARY']  # PLACE FO GRIPS
        incnum = len(odb.steps['Step-1'].frames)
        nel = len(inst.elements) #number of elements
        force = np.zeros((incnum-1, 1))
        i = 0
        j=-1
        Strain = np.zeros((nel*(incnum-1), 3))
        Stress = np.zeros((nel*(incnum-1), 3))
        Elem_intpoint_COORD = np.zeros((nel*(incnum-1), 2))
        EqPlasticStrain = np.zeros((nel*(incnum-1), 1))
        #Extracting RF node labels
        a = nset.nodes
        nodlabel = [a.label - 1 for a in nset.nodes]

        tensi_i = test_config.find('Load_0')
        shear_i = test_config.find('Load_90')

        # Extract data of all time increments
        for frame in odb.steps['Step-1'].frames:
            j=j+1
            if j==0:
                continue
            else:
                if tensi_i == 0:
                    RF = np.array(frame.fieldOutputs['RF'].bulkDataBlocks[0].data.tolist())[nodlabel,0]
                    auxforce = np.sum(RF)
                elif shear_i == 0:
                    RF = np.array(frame.fieldOutputs['RF'].bulkDataBlocks[0].data.tolist())[nodlabel,1]
                    auxforce = np.sum(RF)
                else:
                    RF_x = np.array(frame.fieldOutputs['RF'].bulkDataBlocks[0].data.tolist())[nodlabel,0]
                    RF_y = np.array(frame.fieldOutputs['RF'].bulkDataBlocks[0].data.tolist())[nodlabel,1]
                    auxforce = sqrt((np.sum(RF_x))**2 + (np.sum(RF_y))**2)
                force[j-1, :] = [auxforce]
                for i in range (nel):
                    for k in [x for x in range(len(frame.fieldOutputs['LE'].values[i].data)) if x != 2]:
                        if k==3:
                            Strain[(j-1)*nel+i,k-1] = frame.fieldOutputs['LE'].values[i].data[k] 
                            Stress[(j-1)*nel+i,k-1] = frame.fieldOutputs['S'].values[i].data[k]
                        else:
                            Strain[(j-1)*nel+i,k] = frame.fieldOutputs['LE'].values[i].data[k]
                            Stress[(j-1)*nel+i,k] = frame.fieldOutputs['S'].values[i].data[k]
                    EqPlasticStrain[(j-1)*nel+i,0] = frame.fieldOutputs['SDV1'].values[i].data
                    Elem_intpoint_COORD[(j-1)*nel+i,0] = odb.steps['Step-1'].frames[0].fieldOutputs['COORD'].bulkDataBlocks[0].data[i][0]
                    Elem_intpoint_COORD[(j-1)*nel+i,1] = odb.steps['Step-1'].frames[0].fieldOutputs['COORD'].bulkDataBlocks[0].data[i][1]
                
        data = np.matrix(np.concatenate((Strain, Stress,Elem_intpoint_COORD,EqPlasticStrain), axis=1))
        header = ['Strain_xx','Strain_yy','Strain_xy','Stress_xx','Stress_yy','Stress_xy','IntegrationPoint_CoordX','IntegrationPoint_CoordY','EqPlasticStrain']

        with open(path_and_file_results, 'w') as f:
            writer = csv.writer(f)
            writer = writer.writerow(header)
            for line in data:
                np.savetxt(f, line, fmt='%.6f', delimiter=',   ')

        odb.close()

        with open(path_to_save_force, 'w') as f:
            np.savetxt(f, force, fmt='%.6f')


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
