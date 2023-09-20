"""
    description
    Creation of the input files for MatchID, as well as the extration of the numerical reaction forcefrom Abaqus
        From the .odb file of a 2D test and inp file, it writes the .mesh and .csv files required as input for the synthetic images generation in MatchID

    INPUT: .odb file of the test with the COORD field output request in the nodal position and .inp of the test for element's connectivity
    OUTPUT: .mesh and .csv files with the nodal coordinates of each increment and element's connectivity & NumForce file wth the sum of the reaction force in the grips for each time step

    Run by writting on the command line: "abaqus cae noGUI=AbaqusPythonScriptToMatchID.py"
"""

# IMPORT PACKAGES -------------------------------------------------------------
from PythonScript_main import *
from abaqusConstants import *
from odbAccess import *
import numpy as np
import sys
import os


#  --------------------------------------------------------------- method - MAIN
def main():
    try:
        os.remove('NumForce.dat')
    except:
        pass

    #test_name = 'ArcanSpecimen'  # USE NAME OF TEST
    #test_config = 'Load_90_MatDir_90' #Load_90_MatDir_0 , Load_45_MatDir_0, Load_0_MatDir_0,...

    cwd = os.getcwd()
    cwd_test = os.path.join(cwd,test_config,)
    path_to_save_MatchID = os.path.join(cwd_test,'MatchID', 'FE_DEF')
    path_to_save_force = os.path.join(cwd_test,'MatchID', 'VFM')    
    # ----------open odb with results of the test (it is required the COORD fieldoutput request)
    os.chdir(os.path.join(cwd_test,'FEA'))
    odb_name = test_name + '.odb'
    odb = openOdb(odb_name, readOnly=True)

    errors = odb.diagnosticData.numberOfAnalysisErrors
    if errors != 0:
        return

    # ------------------------------------------------ instance of analysis
    inst = odb.rootAssembly.instances['PART-1-1']
    nset = odb.rootAssembly.instances['PART-1-1'].nodeSets['RIGHTBOUNDARY']  # PLACE FO GRIPS
    incnum = len(odb.steps['Step-1'].frames)
    force = np.zeros((incnum, 1))
    nnod = len(inst.nodes)  # number of nodes
    data = {}
    und_data = {}
    i = 0
    istep = 0
    und_data['CoordZ'] = np.zeros((nnod, 1))
    #CoordZ = 0.0000000000

    #Extracting RF node labels
    a = nset.nodes
    nodlabel = [a.label - 1 for a in nset.nodes]

    #Extracting Surface node labels
    surface = odb.rootAssembly.instances['PART-1-1'].nodeSets['SURFACE'].nodes
    node_surface = [surface.label - 1 for surface in surface]

    #Extracting Surface node labels
    el_surface = odb.rootAssembly.instances['PART-1-1'].elementSets['SURFACE'].elements
    elements_surface = [el_surface.label - 1 for el_surface in el_surface]

    tensi_i = test_config.find('Load_0')
    shear_i = test_config.find('Load_90')

    for frame in odb.steps['Step-1'].frames:  # for every frame/increment
        if i >= istep:
            auxcoord = frame.fieldOutputs['COORD']  # get coordinates data
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
            time = frame.frameValue
            force[i, :] = [auxforce]
            if i == 0:
                und_data['NodeLabel'] = (auxcoord.bulkDataBlocks[0].nodeLabels[node_surface]).astype(int) - 1  # get nodes' label
                und_data['CoordX'] = np.array(auxcoord.bulkDataBlocks[0].data[node_surface])[:, 0]  # get coord 1 for all nodes
                und_data['CoordY'] = np.array(auxcoord.bulkDataBlocks[0].data[node_surface])[:, 1]  # get coord 2 for all nodes
                und_data['CoordZ'] = np.array(auxcoord.bulkDataBlocks[0].data[node_surface])[:, 2]  # get coord 2 for all nodes

                test_K_name = test_name[0:] + '_' + str(i) + '.mesh'
                path_and_file_mesh = os.path.join(path_to_save_MatchID,test_K_name)

                f_0 = open(path_and_file_mesh, 'wb')
                H = ['*Part, name=' + test_name + '\n', '*Node\n']
                f_0.writelines(H)

                for node_number in und_data[
                    'NodeLabel']:  # write undeformed state of test (node label; coordx, coordy, coordz)
                    line = "    {N};    {X};    {Y};    {Z}\n".format(N=und_data['NodeLabel'][node_number] + 1,
                                                                      X=und_data['CoordX'][node_number],
                                                                      Y=und_data['CoordY'][node_number], 
                                                                      Z=und_data['CoordZ'][node_number])
                    f_0.write(line)
                f_0.close()

                El_conect = []
                # Open the file on read only mode
                inp_test_file = test_name + '.inp'
                with open(inp_test_file, 'r') as read_obj:
                    # Read all lines in the file one by one
                    search_begining = '*Element, type='
                    search_end = '*Nset, nset='
                    for line in read_obj:
                        # Check if line contains string
                        if search_begining in line:
                            for line in read_obj:
                                if search_end in line:
                                    break
                                else:
                                    El_conect.append((line.replace(" ", "")).split(','))
                El_conect = [val for j, val in enumerate(El_conect) if j in elements_surface] # Keeping only the elements on the surface element set
                f_0 = open(path_and_file_mesh, 'a')
                H2 = ['*Element\n']
                f_0.writelines(H2)

                for el_conect_line in El_conect:  # write undeformed state of test (node label; coordx, coordy, coordz)
                    line = " {A};  {B};  {C};  {D};  {E}\n".format(A=el_conect_line[0], B=el_conect_line[1],
                                                                 C=el_conect_line[2], D=el_conect_line[3],
                                                                 E=el_conect_line[4])
                    f_0.write(line)
                f_0.close()

            else:
                data['NodeLabel'] = (auxcoord.bulkDataBlocks[0].nodeLabels[node_surface]).astype(int) - 1  # get nodes' label
                data['CoordX'] = np.array(auxcoord.bulkDataBlocks[0].data[node_surface])[:, 0]  # get coord 1 for all nodes
                data['CoordY'] = np.array(auxcoord.bulkDataBlocks[0].data[node_surface])[:, 1]  # get coord 2 for all nodes
                data['CoordZ'] = np.array(auxcoord.bulkDataBlocks[0].data[node_surface])[:, 2]  # get coord 2 for all nodes

                test_K_name = test_name[0:] + '_' + str(i) + '.csv'
                path_and_file_csv = os.path.join(path_to_save_MatchID, test_K_name)

                f = open(path_and_file_csv, 'wb')
                f.writelines(H)

                for node_number in data['NodeLabel']:  # write deformed state of test (node label; coordx, coordy, coordz)
                    line = "    {N};    {X};    {Y};    {Z}\n".format(N=data['NodeLabel'][node_number] + 1,
                                                                      X=data['CoordX'][node_number],
                                                                      Y=data['CoordY'][node_number], 
                                                                      Z=data['CoordZ'][node_number])
                    f.write(line)
                f.close()
        i += 1

    odb.close()
    
    path_and_file_force = os.path.join(path_to_save_force,'NumForce.csv')

    with open(path_and_file_force, 'w') as f:
        np.savetxt(f, np.delete(force,0,0), fmt='%.6f')



# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
