"""
    description
   Extraction of the Major and minor strain components

    INPUT: .odb file of the test with the LE, SDV1, S and EVOL field output request 
    OUTPUT: .csv files with principal stresses and strains, eq. plastic strain, and element volume


    Run by writting on the command line: "abaqus cae noGUI=AbaqusPythonScript_Results.py"
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

    #test_name = 'ArcanSpecimen'  # USE NAME OF TEST
    #test_config = 'Load_90_MatDir_90' #Load_90_MatDir_0 , Load_45_MatDir_0, Load_0_MatDir_0,...

    cwd = os.getcwd()
    cwd_test = os.path.join(cwd,test_config)
    #path_to_save_Results = os.path.join(cwd,'\\Results\\FEA_results\\')
    #path_to_save_Results = 'D:/JDHenriques/ArcanTest/ROI/Load_90_MatDir_0/FEA/Results'
    # ----------open odb with results of the test (it is required the COORD fieldoutput request)
    os.chdir(os.path.join(cwd_test,'FEA'))
    odb_name = test_name + '.odb'
    odb = openOdb(odb_name, readOnly=True)

    errors = odb.diagnosticData.numberOfAnalysisErrors
    if errors != 0:
        return

    # ------------------------------------------------ instance of analysis
    inst = odb.rootAssembly.instances['PART-1-1']
    incnum = len(odb.steps['Step-1'].frames)
    nel = len(inst.elementSets['SURFACE'].elements) #number of elements
    i = 0
    minPrincipalStrain = np.zeros((nel, 1))
    maxPrincipalStrain = np.zeros((nel, 1))
    eqplasticstrain = np.zeros((nel, 1))
    minPrincipalStress = np.zeros((nel, 1))
    maxPrincipalStress = np.zeros((nel, 1))
    ElemVolume = np.zeros((nel, 1))
    Elem_intpoint_COORD = np.zeros((nel, 2))
    Stress = np.zeros((nel, 3))
    PressureStress = np.zeros((nel, 1))
    MisesStress = np.zeros((nel, 1))
    ThirdInvStress = np.zeros((nel, 1))

    #Retrieving the connectivity of the elements to calculate the coordinates of the centroid
    element_conec = [inst.elementSets['SURFACE'].elements[j].connectivity for j in range(0,nel)]
    
    # Pre-calculate fieldOutputs for Step-1 frames[-1] to avoid redundant calls
    field_outputs_LE = odb.steps['Step-1'].frames[-1].fieldOutputs['LE'].getSubset(position=CENTROID, region=inst.elementSets['SURFACE']).values
    field_outputs_SDV1 = odb.steps['Step-1'].frames[-1].fieldOutputs['SDV1'].getSubset(position=CENTROID, region=inst.elementSets['SURFACE']).values
    field_outputs_S = odb.steps['Step-1'].frames[-1].fieldOutputs['S'].getSubset(position=CENTROID, region=inst.elementSets['SURFACE']).values
    field_outputs_EVOL = odb.steps['Step-1'].frames[-1].fieldOutputs['EVOL'].values
    field_outputs_COORD = odb.steps['Step-1'].frames[0].fieldOutputs['COORD'].bulkDataBlocks[0].data

    # loop 
    for i in range(nel):
        # Get values for 'LE', 'SDV1', 'S', and 'EVOL' only once for each element
        le_values = field_outputs_LE[i]
        s_values = field_outputs_S[i]
        
        ElemVolume[i] = field_outputs_EVOL[i].data
        eqplasticstrain[i] = field_outputs_SDV1[i].data
        
        minPrincipalStrain[i] = le_values.minPrincipal
        maxPrincipalStrain[i] = le_values.maxPrincipal

        minPrincipalStress[i] = s_values.minPrincipal
        maxPrincipalStress[i] = s_values.maxPrincipal
        
        sum_x = 0
        sum_y = 0
        for k in range(len(element_conec[0])):
            elem_index = element_conec[i][k] - 1
            sum_x += field_outputs_COORD[elem_index][0]
            sum_y += field_outputs_COORD[elem_index][1]

        Elem_intpoint_COORD[i, 0] = sum_x / len(element_conec[0])
        Elem_intpoint_COORD[i, 1] = sum_y / len(element_conec[0])

        Stress[i, 0] = s_values.data[0]
        Stress[i, 1] = s_values.data[1]
        Stress[i, 2] = s_values.data[3]
        PressureStress[i] = s_values.press
        MisesStress[i] = s_values.mises
        ThirdInvStress[i] = s_values.inv3

    StressTriax = -(PressureStress/MisesStress)

    data = np.matrix(np.concatenate((minPrincipalStrain, maxPrincipalStrain, eqplasticstrain, minPrincipalStress, maxPrincipalStress, ElemVolume, Elem_intpoint_COORD, Stress, StressTriax, MisesStress, ThirdInvStress), axis=1))
    header = ['MinPrincipalStrain','MaxPrincipalStrain','EqPlasticStrain','MinPrincipalStress','MaxPrincipalStress','ElemVolume','IntegrationPoint_CoordX','IntegrationPoint_CoordY','Stress_xx','Stress_yy','Stress_xy', 'StressTriaxility', 'MisesStress', 'ThirdInvStress']
    
    path_and_file_results = (cwd + '\\Results\\FEA_results\\Strain_Stress_' + test_config + '.csv')

    with open(path_and_file_results, 'w') as f:
       writer = csv.writer(f)
       writer = writer.writerow(header)
       for line in data:
           np.savetxt(f, line, fmt='%.6f', delimiter=',   ')

    #Extract results for all increments 
    i = 0
    j= -1
    minPrincipalStress = np.zeros((nel*(incnum-1), 1))
    maxPrincipalStress = np.zeros((nel*(incnum-1), 1))
    minPrincipalStrain = np.zeros((nel*(incnum-1), 1))
    maxPrincipalStrain = np.zeros((nel*(incnum-1), 1))
    Strain_ref = np.zeros((nel*(incnum-1), 3))
    Stress_ref = np.zeros((nel*(incnum-1), 3))
    EqPlasticStrain_ref = np.zeros((nel*(incnum-1), 1))
    
    excluded_indices = {2, 4, 5}  # Set of indices to exclude
    
    for frame in odb.steps['Step-1'].frames:
        j = j+1
        if j==0:
            continue
        else:
            element_LE = frame.fieldOutputs['LE'].getSubset(position=CENTROID, region=inst.elementSets['SURFACE']).values
            element_S = frame.fieldOutputs['S'].getSubset(position=CENTROID, region=inst.elementSets['SURFACE']).values
            element_SDV1 = frame.fieldOutputs['SDV1'].getSubset(position=CENTROID, region=inst.elementSets['SURFACE']).values
        
            for i in range(nel):
                # Update minPrincipalStress and maxPrincipalStress arrays
                minPrincipalStress[(j - 1) * nel + i] = element_S[i].minPrincipal
                maxPrincipalStress[(j - 1) * nel + i] = element_S[i].maxPrincipal
            
                # Update minPrincipalStrain and maxPrincipalStrain arrays
                minPrincipalStrain[(j - 1) * nel + i] = element_LE[i].minPrincipal
                maxPrincipalStrain[(j - 1) * nel + i] = element_LE[i].maxPrincipal
            
                # Update EqPlasticStrain_ref array (assuming it's a scalar value)
                EqPlasticStrain_ref[(j - 1) * nel + i, 0] = element_SDV1[i].data
        
                # Update Strain_ref and Stress_ref arrays
                for k in [x for x in range(len(frame.fieldOutputs['LE'].values[i].data)) if x not in excluded_indices]:
                    if k==3:
                        Strain_ref[(j-1)*nel+i,k-1] = element_LE[i].data[k]
                        Stress_ref[(j-1)*nel+i,k-1] = element_S[i].data[k]
                    else:
                        Strain_ref[(j-1)*nel+i,k] = element_LE[i].data[k]
                        Stress_ref[(j-1)*nel+i,k] = element_S[i].data[k]
                      
    data_ref = np.matrix(np.concatenate((Strain_ref, Stress_ref,EqPlasticStrain_ref,minPrincipalStrain,maxPrincipalStrain,minPrincipalStress,maxPrincipalStress), axis=1))
    header_ref = ['Strain_xx','Strain_yy','Strain_xy','Stress_xx','Stress_yy','Stress_xy','EqPlasticStrain','MinPrincipalStrain','MaxPrincipalStrain','MinPrincipalStress','MaxPrincipaStress']

    path_and_file_results_ref = (cwd + '\\Results\\FEA_results\\' + test_config + '_AllSteps.csv')
    with open(path_and_file_results_ref, 'w') as f:
        writer = csv.writer(f)
        writer = writer.writerow(header_ref)
        for line in data_ref:
            np.savetxt(f, line, fmt='%.6f', delimiter=',   ')

    odb.close()

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
