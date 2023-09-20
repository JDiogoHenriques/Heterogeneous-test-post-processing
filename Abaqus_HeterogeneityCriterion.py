"""
    description
   Plots principal stress and strain components

    INPUT: .csv file of the test with the LE, SDV1, S and EVOL field output request 
    OUTPUT: Plots the principal stress and strain components


    Run by writting on the command line: "python Abaqus_HeterogeneityCriterion.py"
"""

# IMPORT PACKAGES -------------------------------------------------------------
from PythonScript_main import *
import numpy as np
import os
import csv
import matplotlib.pyplot as plt

#  --------------------------------------------------------------- method - MAIN
def main():

    #test_name = 'ArcanSpecimen'  # USE NAME OF TEST
    #test_config = 'Load_0_MatDir_0' #Load_90_MatDir_0 , Load_45_MatDir_0, Load_0_MatDir_0,...

    # Heterogeinity Criterion weights
    wa1 = 1.00
    wa2 = 4.00
    wa3 = 0.25
    wa4 = 0.80
    wa5 = 0.40
    wr1 = 0.13
    wr2 = 0.02
    wr3 = 0.25
    wr4 = 0.35
    wr5 = 0.25
    #Strain states interval 
    Interval_PrincStrains = 0.03

    cwd = os.getcwd()
    cwd_test = os.path.join(cwd,test_config,)
    path_to_read_Results = os.path.join(cwd,'Results\\FEA_results\\')
    res_file_name = 'Strain_Stress_'

    #Importing and data management
    all_results = []

    files = [filename for filename in os.listdir(path_to_read_Results) if
             filename.startswith(res_file_name) and filename.endswith(".csv")]
    if not files:
        print('ERROR: No .csv files with FEA results.')
        exit()
    
    i = 0
    for file_name in files:  # read each file refering to a test
        with open(path_to_read_Results + file_name, 'r', errors='ignore') as infile:
            reader = list(csv.reader(infile, delimiter=','))
            all_results.append(np.array(reader[2:]))
            if i==0:
                n_points=np.array([len(all_results[i])])  # number of data points of first test
            else:
                n_points = np.append(n_points,len(all_results[i]))  # number of data points for each test
            i = i + 1  # counter

    #--------------- Heterogeneity criterion computation ------------------------------   
    MinPrincStrain = np.zeros((i,n_points[0],1))
    MaxPrincStrain = np.zeros((i,n_points[0],1))
    EqPlastStrain = np.zeros((i,n_points[0],1))
    ElemVolume = np.zeros((i,n_points[0],1))
    Std_PrincStrains = np.zeros((i,1,1))
    Range_PrincStrains = np.zeros((i,1,1))
    Std_EqPlastStrain = np.zeros((i,1,1))
    Av_EqPlastStrain = np.zeros((i,1,1))
    Max_EqPlasStrain_uni_tens = np.zeros((i,1,1))
    Max_EqPlasStrain_uni_shear = np.zeros((i,1,1))
    Max_EqPlasStrain_uni_planestraintension = np.zeros((i,1,1))
    Max_EqPlasStrain_uni_equibiaxial = np.zeros((i,1,1))
    Max_EqPlasStrain_uni_compression = np.zeros((i,1,1)) 
    Max_test_EqPlasticStrain = np.zeros((i,1,1)) 
    Max_EqPlasticStrain = np.zeros((i,1,1)) 
    IT1_aux = np.zeros((i,1,1)) 
    IT1 = np.zeros((i,1)) 

    np.seterr(invalid='ignore')

    for test in range(0, i):
        for line in range(0,n_points[0]):
            MinPrincStrain[test,line,0] = all_results[test][line][0]
            MaxPrincStrain[test,line,0] = all_results[test][line][1]
            EqPlastStrain[test,line,0] = all_results[test][line][2]
            ElemVolume[test,line,0] = all_results[test][line][5]   
        Std_PrincStrains[test,0,0] = np.std((MinPrincStrain/MaxPrincStrain)[test])
        Range_PrincStrains[test,0,0] = max((MinPrincStrain/MaxPrincStrain)[test])-min((MinPrincStrain/MaxPrincStrain)[test])
        Std_EqPlastStrain[test,0,0] = np.std((EqPlastStrain)[test])
        Av_EqPlastStrain[test,0,0] = sum((EqPlastStrain*ElemVolume)[test])/sum((ElemVolume)[test])

    # Strain states calculation
    States_PrincStrain = MinPrincStrain/MaxPrincStrain

    # Strain uniaxial tension
    Range_uni_tens = [-0.5-Interval_PrincStrains, -0.5+Interval_PrincStrains]
    index_max_uni_tens = np.where((States_PrincStrain < Range_uni_tens[1]) & (States_PrincStrain > Range_uni_tens[0]))[1]
    tens_test = np.where((States_PrincStrain < Range_uni_tens[1]) & (States_PrincStrain > Range_uni_tens[0]))[0]

    for index in range(0, len(available_test_configs)**2):
        ind1 = index_max_uni_tens[np.where(tens_test==index)[0]]      
        if ind1.size>=1:
            Max_EqPlasStrain_uni_tens[index,0,0] = max(abs(EqPlastStrain[index,ind1,0]))
        else:
            Max_EqPlasStrain_uni_tens[index,0,0] = 0

    #Strain pure shear
    Range_pure_shear = [-1-Interval_PrincStrains, -1+Interval_PrincStrains]
    index_max_uni_shear = np.where((States_PrincStrain < Range_pure_shear[1]) & (States_PrincStrain > Range_pure_shear[0]))[1]
    shear_test = np.where((States_PrincStrain < Range_pure_shear[1]) & (States_PrincStrain > Range_pure_shear[0]))[0]

    for index1 in range(0, len(available_test_configs)**2):
        ind2 = index_max_uni_shear[np.where(shear_test==index1)[0]]      
        if ind2.size>=1:
            Max_EqPlasStrain_uni_shear[index1,0,0] = max(abs(EqPlastStrain[index1,ind2,0]))
        else:
            Max_EqPlasStrain_uni_shear[index1,0,0] = 0

    #Strain planestrain tension
    Range_planestraintension = [0-Interval_PrincStrains, 0+Interval_PrincStrains]
    index_max_planestraintension = np.where((States_PrincStrain < Range_planestraintension[1]) & (States_PrincStrain > Range_planestraintension[0]))[1]
    planestrain_test = np.where((States_PrincStrain < Range_planestraintension[1]) & (States_PrincStrain > Range_planestraintension[0]))[0]

    for index2 in range(0, len(available_test_configs)**2):
        ind3 = index_max_planestraintension[np.where(planestrain_test==index2)[0]]      
        if ind3.size>=1:
            Max_EqPlasStrain_uni_planestraintension[index2,0,0] = max(abs(EqPlastStrain[index2,ind3,0]))
        else:
            Max_EqPlasStrain_uni_planestraintension[index2,0,0] = 0        

    #Strain equibiaxial tension
    Range_equibiaxial = [1-Interval_PrincStrains, 1+Interval_PrincStrains]
    index_max_equibiaxial = np.where((States_PrincStrain < Range_equibiaxial[1]) & (States_PrincStrain > Range_equibiaxial[0]))[1]
    equibiaxial_test = np.where((States_PrincStrain < Range_equibiaxial[1]) & (States_PrincStrain > Range_equibiaxial[0]))[0]

    for index3 in range(0, len(available_test_configs)**2):
        ind4 = index_max_equibiaxial[np.where(equibiaxial_test==index3)[0]]      
        if ind4.size>=1:
            Max_EqPlasStrain_uni_equibiaxial[index3,0,0] = max(abs(EqPlastStrain[index3,ind4,0]))
        else:
            Max_EqPlasStrain_uni_equibiaxial[index3,0,0] = 0   

    #Strain Uniaxial compression
    Range_uni_compression = [-2-Interval_PrincStrains, -2+Interval_PrincStrains]
    index_max_compression = np.where((States_PrincStrain < Range_uni_compression[1]) & (States_PrincStrain > Range_uni_compression[0]))[1]
    compression_test = np.where((States_PrincStrain < Range_uni_compression[1]) & (States_PrincStrain > Range_uni_compression[0]))[0]

    for index4 in range(0, len(available_test_configs)**2):
        ind5 = index_max_compression[np.where(compression_test==index4)[0]]      
        if ind5.size>=1:
            Max_EqPlasStrain_uni_compression[index4,0,0] = max(abs(EqPlastStrain[index4,ind5,0]))
        else:
            Max_EqPlasStrain_uni_compression[index4,0,0] = 0    

    # Max from Test and heterogeneity criterion computation calculation 
    for index5 in range(0, len(available_test_configs)**2):
        Max_test_EqPlasticStrain[index5,0,0] = max(abs(EqPlastStrain[index5]))
        Max_EqPlasticStrain[index5,0,0] = (Max_test_EqPlasticStrain[index5] + Max_EqPlasStrain_uni_compression[index5] + \
            Max_EqPlasStrain_uni_equibiaxial[index5] + Max_EqPlasStrain_uni_planestraintension[index5] + \
            Max_EqPlasStrain_uni_shear[index5] + Max_EqPlasStrain_uni_tens[index5])/6
        IT1_aux[index5,0,0] = wr1*(Std_PrincStrains[index5]/wa1) + wr2*(Range_PrincStrains[index5]/wa2) + \
            wr3*(Std_EqPlastStrain[index5]/wa3) + wr4*(Max_EqPlasticStrain[index5]/wa4) + wr5*(Av_EqPlastStrain[index5]/wa5)
        IT1[index5,0] = IT1_aux [index5]

    #Plot data preparation
    IT1 = np.reshape(IT1,(len(available_test_configs),len(available_test_configs)), order='F')
    Range_PrincStrains_grid = np.reshape(Range_PrincStrains,(len(available_test_configs),len(available_test_configs)), order='F')
    Std_PrincStrains_grid = np.reshape(Std_PrincStrains,(len(available_test_configs),len(available_test_configs)), order='F')
    Std_EqPlastStrain_grid = np.reshape(Std_EqPlastStrain,(len(available_test_configs),len(available_test_configs)), order='F') 
    Max_EqPlasticStrain_grid = np.reshape(Max_EqPlasticStrain,(len(available_test_configs),len(available_test_configs)), order='F') 
    Av_EqPlastStrain_grid = np.reshape(Av_EqPlastStrain,(len(available_test_configs),len(available_test_configs)), order='F') 
    test_configs = np.array([float(x) for x in available_test_configs]) 
    Load_angle = test_configs
    Material_Direction = test_configs
    Load_angle_arr = np.repeat(test_configs,len(available_test_configs))
    Material_Direction_arr = np.transpose(np.tile(test_configs,(1,len(available_test_configs))))
    Load_angle_grid = np.reshape(Load_angle_arr,(len(available_test_configs),len(available_test_configs)), order='F')
    Material_Direction_grid = np.reshape(Material_Direction_arr,(len(available_test_configs),len(available_test_configs)), order='F')
    # ------------- graph properties ------------- 
    FS = 24
    plt.rcParams['axes.facecolor'] = (1, 1, 1)
    plt.rcParams['figure.facecolor'] = (1, 1, 1)
    plt.rcParams["font.family"] = "sans"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams['font.size'] = FS
    params = {"ytick.color" : (0, 0, 0),
          "xtick.color" : (0, 0, 0),
          "grid.color" : (0, 0 , 0),
          "text.color" : (0, 0, 0),
          "axes.labelcolor" : (0, 0, 0),
          "axes.edgecolor" : (.15, .15, .15),
          "text.usetex": False}
    plt.rcParams.update(params)

    # -------------Heterogeneity criterion plot------------- 
    fig1 = plt.figure(figsize=(18,9))
    ax = plt.axes([0.12, 0.14, .86, .83])

    extents = (-90/(len(available_test_configs)+1),90*(1+(1/(len(available_test_configs)+1))),\
        -90/(len(available_test_configs)+1),90*(1+(1/(len(available_test_configs)+1))))

    plt.imshow(IT1,cmap='Reds',origin='lower',extent=extents)
    plt.colorbar()
    plt.title("Heterogeinity criterion")
    plt.xlabel("Load angle, "+r'$\degree$', fontsize=FS, labelpad=10)
    plt.ylabel("Material direction, "+r'$\degree$', fontsize=FS, labelpad=10)
    ax.set_xticks(test_configs)
    ax.set_yticks(test_configs)

    fig1.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'HeterCriterion_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    # -------------Heterogeneity criterion 2D plot------------- 
    fig2 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Load_angle,IT1[0,0:], color='red', marker='v', zorder=10, markersize=12, clip_on=False, label=r'$\theta=0\degree$', linewidth=3) # Mat Dir 0
    plt.plot(Load_angle,IT1[1,0:], color='black', marker='o', zorder=10, markersize=12, clip_on=False, label=r'$\theta=45\degree$', linewidth=3) # Mat Dir 45
    plt.plot(Load_angle,IT1[2,0:], color='blue', marker='^', zorder=10, markersize=12, clip_on=False, label=r'$\theta=90\degree$', linewidth=3) # Mat Dir 90
    
    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Heterogeinity criterion")
    plt.xlabel(r'$\alpha\ [\degree]$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$IT_{1}$', fontsize=FS, labelpad=10)
    plt.xlim(-0.5, 90.5)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    ax.set_xticks(Load_angle)

    fig2.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'HeterCriterion_2Dplot_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    # -------------Strain range 2D plot------------- 
    fig3 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Load_angle,Range_PrincStrains_grid[0,0:], color='red', marker='v', zorder=10, markersize=12, clip_on=False, label=r'$\theta=0\degree$', linewidth=3) # Mat Dir 0
    plt.plot(Load_angle,Range_PrincStrains_grid[1,0:], color='black', marker='o', zorder=10, markersize=12, clip_on=False, label=r'$\theta=45\degree$', linewidth=3) # Mat Dir 45
    plt.plot(Load_angle,Range_PrincStrains_grid[2,0:], color='blue', marker='^', zorder=10, markersize=12, clip_on=False, label=r'$\theta=90\degree$', linewidth=3) # Mat Dir 90
    
    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Principal strain range")
    plt.xlabel(r'$\alpha\ [\degree]$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$(\epsilon_2/\epsilon_1)_\mathrm{R}$', fontsize=FS, labelpad=10)
    plt.xlim(-0.5, 90.5)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    ax.set_xticks(Load_angle)

    fig3.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'PrincStrainRange_2Dplot_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    # -------------Strain std 2D plot------------- 
    fig3 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Load_angle,Std_PrincStrains_grid[0,0:], color='red', marker='v', zorder=10, markersize=12, clip_on=False, label=r'$\theta=0\degree$', linewidth=3) # Mat Dir 0
    plt.plot(Load_angle,Std_PrincStrains_grid[1,0:], color='black', marker='o', zorder=10, markersize=12, clip_on=False, label=r'$\theta=45\degree$', linewidth=3) # Mat Dir 45
    plt.plot(Load_angle,Std_PrincStrains_grid[2,0:], color='blue', marker='^', zorder=10, markersize=12, clip_on=False, label=r'$\theta=90\degree$', linewidth=3) # Mat Dir 90
    
    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Principal strain standard deviation")
    plt.xlabel(r'$\alpha\ [\degree]$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\mathrm{Std}(\epsilon_2/\epsilon_1)$', fontsize=FS, labelpad=10)
    plt.xlim(-0.5, 90.5)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    ax.set_xticks(Load_angle)

    fig3.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'PrincStrainStd_2Dplot_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    # -------------Plastic strain std 2D plot------------- 
    fig4 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Load_angle,Std_EqPlastStrain_grid[0,0:], color='red', marker='v', zorder=10, markersize=12, clip_on=False, label=r'$\theta=0\degree$', linewidth=3) # Mat Dir 0
    plt.plot(Load_angle,Std_EqPlastStrain_grid[1,0:], color='black', marker='o',zorder=10, markersize=12, clip_on=False, label=r'$\theta=45\degree$', linewidth=3) # Mat Dir 45
    plt.plot(Load_angle,Std_EqPlastStrain_grid[2,0:], color='blue', marker='^',zorder=10, markersize=12, clip_on=False, label=r'$\theta=90\degree$', linewidth=3) # Mat Dir 90
    
    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Plastic strain standard deviation")
    plt.xlabel(r'$\alpha\ [\degree]$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\mathrm{Std}(\bar{\epsilon}^{\mathrm{P}})$', fontsize=FS, labelpad=10)
    plt.xlim(-0.5, 90.5)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    ax.set_xticks(Load_angle)

    fig4.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'PlasticStrainStd_2Dplot_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    # -------------Plastic strain max 2D plot------------- 
    fig5 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Load_angle,Max_EqPlasticStrain_grid[0,0:], color='red', marker='v', zorder=10, markersize=12, clip_on=False, label=r'$\theta=0\degree$', linewidth=3) # Mat Dir 0
    plt.plot(Load_angle,Max_EqPlasticStrain_grid[1,0:], color='black', marker='o', zorder=10, markersize=12, clip_on=False, label=r'$\theta=45\degree$', linewidth=3) # Mat Dir 45
    plt.plot(Load_angle,Max_EqPlasticStrain_grid[2,0:], color='blue', marker='^', zorder=10, markersize=12, clip_on=False, label=r'$\theta=90\degree$', linewidth=3) # Mat Dir 90
    
    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Plastic strain max value")
    plt.xlabel(r'$\alpha\ [\degree]$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\bar{\epsilon}^{\mathrm{P}}_{\mathrm{Max}}$', fontsize=FS, labelpad=10)
    plt.xlim(-0.5, 90.5)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    ax.set_xticks(Load_angle)

    fig5.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'PlasticStrainMax_2Dplot_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    # -------------Plastic strain average 2D plot------------- 
    fig6 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Load_angle,Av_EqPlastStrain_grid[0,0:], color='red', marker='v', zorder=10, markersize=12, clip_on=False, label=r'$\theta=0\degree$', linewidth=3) # Mat Dir 0
    plt.plot(Load_angle,Av_EqPlastStrain_grid[1,0:], color='black', marker='o', zorder=10, markersize=12, clip_on=False, label=r'$\theta=45\degree$', linewidth=3) # Mat Dir 45
    plt.plot(Load_angle,Av_EqPlastStrain_grid[2,0:], color='blue', marker='^', zorder=10, markersize=12, clip_on=False, label=r'$\theta=90\degree$', linewidth=3) # Mat Dir 90
    
    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Plastic strain average value")
    plt.xlabel(r'$\alpha\ [\degree]$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\mathrm{Av}(\bar{\epsilon}^{\mathrm{P}})$', fontsize=FS, labelpad=10)
    plt.xlim(-0.5, 90.5)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    ax.set_xticks(Load_angle)

    fig6.savefig(cwd + '\\Results\\FEA_heter_criterion\\' + 'PlasticStrainAverage_2Dplot_' + str(len(available_test_configs)) + 'tests.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    if show_plots == 'yes':
        plt.show()
    else:
        pass
    
    index_heter = np.unravel_index(np.argmax(IT1, axis=None), IT1.shape)

    print("The most heterogeneous configuration according to the criterion is: Load="+str(int(Load_angle[2]))+" MatDir="+str(int(Material_Direction[2])))
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()