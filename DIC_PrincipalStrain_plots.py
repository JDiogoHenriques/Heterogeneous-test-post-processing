"""
    description
   Plots principal stress and strain components

    INPUT: .csv file of the test with the LE, SDV1, S and EVOL field output request 
    OUTPUT: Plots the principal stress and strain components


    Run by writting on the command line: "python DIC_Principaltrain_plots.py"
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
    cwd = os.getcwd()
    cwd_test = os.path.join(cwd,test_config,)
    path_to_save_Results = os.path.join(cwd,'Results\\FEA_results\\')
    path_DIC = os.path.join(cwd_test,'MatchID\\DIC\\')

    #Importing and data management
    all_results = []
    with open(path_to_save_Results + 'Strain_Stress_'+ test_config + '.csv', 'r') as csv_file:
        reader = list(csv.reader(csv_file, delimiter=','))
        #all_results.append({header: list(map(float, values)) for header, *values in zip(*reader)})
        n_points = len(reader)
        all_results = np.array(reader[2:n_points])

    MinPrincStrain = all_results[:,0]
    MaxPrincStrain = all_results[:,1]
    EqPlastStrain = all_results[:,2]
    MinPrincStress = all_results[:,3]
    MaxPrincStress = all_results[:,4]

    MinPrincStrain = np.array([float(x) for x in MinPrincStrain])
    MaxPrincStrain = np.array([float(x) for x in MaxPrincStrain])
    EqPlastStrain = np.array([float(x) for x in EqPlastStrain])
    MinPrincStress = np.array([float(x) for x in MinPrincStress])
    MaxPrincStress = np.array([float(x) for x in MaxPrincStress])
    
    ElasticPoints = np.array([x for x in range(len(EqPlastStrain)) if EqPlastStrain[x]==0])

    if ElasticPoints.size>0:
        ElastMinPrincStrain = np.array(MinPrincStrain[ElasticPoints])
        ElastMaxPrincStrain = np.array(MaxPrincStrain[ElasticPoints])
        ElastMinPrincStress = np.array(MinPrincStress[ElasticPoints])
        ElastMaxPrincStress = np.array(MaxPrincStress[ElasticPoints])

        PlastMinPrincStrain = np.delete(MinPrincStrain, ElasticPoints, 0)
        PlastMaxPrincStrain = np.delete(MaxPrincStrain, ElasticPoints, 0)
        PlastMinPrincStress = np.delete(MinPrincStress, ElasticPoints, 0)
        PlastMaxPrincStress = np.delete(MaxPrincStress, ElasticPoints, 0)
    else:
        pass 

    #Importing and data management
    DIC_results = []
    with open(path_DIC + 'Image_0000_0_Numerical_40_0.def' + '.csv', 'r') as csv_file:
        reader_DIC = list(csv.reader(csv_file, delimiter=';'))
        #all_results.append({header: list(map(float, values)) for header, *values in zip(*reader)})
        n_points_DIC = len(reader_DIC)
        DIC_results = np.array(reader_DIC[1:n_points_DIC])   

    MinPrincStrain_DIC = DIC_results[:,10]
    MaxPrincStrain_DIC = DIC_results[:,9]

    MinPrincStrain_DIC = np.array([float(x) for x in MinPrincStrain_DIC])
    MaxPrincStrain_DIC = np.array([float(x) for x in MaxPrincStrain_DIC])

#--------------- Strain processing ------------------------------
#Strain Uniaxial compression
    Strain_uni_comp_x = np.linspace(-10,0)
    Strain_uni_comp_y = Strain_uni_comp_x/-2

#Strain pure shear
    Strain_pure_shear_x = np.linspace(-10,0)
    Strain_pure_shear_y = Strain_pure_shear_x/-1

#Strain uniaxial tension
    Strain_uni_tens_x = np.linspace(-10,0)
    Strain_uni_tens_y = Strain_uni_tens_x/-0.5

#Strain planestrain tension
    Strain_planestra_tens_y = np.linspace(0,10)
    Strain_planestra_tens_x = np.zeros(len(Strain_planestra_tens_y))

#Strain equibiaxial tension
    Strain_equibiaxial_tens_x = np.linspace(0,10)
    Strain_equibiaxial_tens_y = Strain_equibiaxial_tens_x

#Forming Limit Diagram
    FLC_Eps2_neg = np.linspace(-10,0)
    FLC_Eps2_pos = np.linspace(0.000000001,10)
    FLC_Eps1_neg = -FLC_Eps2_neg+0.194
    FLC_Eps1_pos = 0.5*FLC_Eps2_pos+0.194
    FLC_Eps1 = np.column_stack((FLC_Eps1_neg,FLC_Eps1_pos))
    FLC_Eps2 = np.column_stack((FLC_Eps2_neg,FLC_Eps2_pos))

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

# -------------Principal strains graph------------- 
    fig1 = plt.figure(figsize=(18,9))
    ax = plt.axes([0.12, 0.14, .86, .83])

    colors = EqPlastStrain
    if ElasticPoints.size>0:
        colors = np.delete(colors, ElasticPoints, 0)
        elastic = plt.scatter(ElastMinPrincStrain,ElastMaxPrincStrain, s=35, c = 'grey', marker = '+', label=r'Elastic')
        plastic = plt.scatter(PlastMinPrincStrain,PlastMaxPrincStrain, c = colors, cmap='jet', s=15, label=r'Plastic')
    else:
        elastic = plt.scatter([],[], s=35, c = 'grey', marker = '+', label=r'Elastic')
        plastic = plt.scatter(MinPrincStrain,MaxPrincStrain, c = colors, cmap='jet', s=15, label=r'Plastic')
  
    color_bar= plt.colorbar()

    DIC = plt.scatter(MinPrincStrain_DIC,MaxPrincStrain_DIC, c = 'grey', marker = 'o',  s=10, alpha=0.2, label=r'DIC')

    plt.plot(Strain_uni_comp_x, Strain_uni_comp_y, color='black', linewidth=0.8)
    plt.plot(Strain_pure_shear_x, Strain_pure_shear_y, color='black', linewidth=0.8)
    plt.plot(Strain_uni_tens_x,Strain_uni_tens_y, color='black', linewidth=0.8)
    plt.plot(Strain_planestra_tens_x,Strain_planestra_tens_y, color='black', linewidth=0.8)
    plt.plot(Strain_equibiaxial_tens_x,Strain_equibiaxial_tens_y, color='black', linewidth=0.8)
    plt.plot(FLC_Eps2,FLC_Eps1, color='grey', linewidth=0.7)

    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1,
                handleheight=1.5, labelspacing=0.5)
    plt.title("Test configuration: " + test_config)
    plt.xlabel("Minor principal strain, "+r'$\epsilon_2$', fontsize=FS, labelpad=10)
    plt.ylabel("Major principal strain, "+r'$\epsilon_1$', fontsize=FS, labelpad=10)
    plt.xlim(-max(-abs(min(MinPrincStrain))*1.05, 0.37), max(abs(min(MinPrincStrain))*1.05, 0.37))
    plt.ylim(0, max(abs(max(MaxPrincStrain))*1.05, 0.39))
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
  


    ax.annotate('Uniaxial \ncompression',xy=(0.145, 0.475), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Pure \nshear',xy=(0.13, 0.86), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Uniaxial \ntension',xy=(0.26, 0.868), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Plane strain \ntension',xy=(0.453, 0.87), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Equibiaxial \ntension',xy=(0.71, 0.7), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('FLC',xy=(0.568, 0.65), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))        
    ax.annotate(r'$\bar{\epsilon}^\mathrm{p}$',xy=(0.855, 0.96), xycoords='figure fraction', fontsize=28, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))    
    ax.annotate('b',xy=(0.02, 0.96), xycoords='figure fraction', fontsize=28, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))    

    fig1.savefig(cwd + '\\Results\\FEA_stress_strain_plots\\' + test_config + '\\' + 'DIC_FEA_Compar_PrincStrain_' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

# -------------Principal strains graph------------- 
    fig2 = plt.figure(figsize=(18,9))
    ax = plt.axes([0.12, 0.14, .86, .83])

    DIC = plt.scatter(MinPrincStrain_DIC,MaxPrincStrain_DIC, c = 'grey', marker = 'o',  s=10, alpha=0.4, label=r'DIC')

    plt.plot(Strain_uni_comp_x, Strain_uni_comp_y, color='black', linewidth=0.8)
    plt.plot(Strain_pure_shear_x, Strain_pure_shear_y, color='black', linewidth=0.8)
    plt.plot(Strain_uni_tens_x,Strain_uni_tens_y, color='black', linewidth=0.8)
    plt.plot(Strain_planestra_tens_x,Strain_planestra_tens_y, color='black', linewidth=0.8)
    plt.plot(Strain_equibiaxial_tens_x,Strain_equibiaxial_tens_y, color='black', linewidth=0.8)
    plt.plot(FLC_Eps2,FLC_Eps1, color='grey', linewidth=0.7)

    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1,
                handleheight=1.5, labelspacing=0.5)
    plt.title("Test configuration: " + test_config)
    plt.xlabel("Minor principal strain, "+r'$\epsilon_2$', fontsize=FS, labelpad=10)
    plt.ylabel("Major principal strain, "+r'$\epsilon_1$', fontsize=FS, labelpad=10)
    plt.xlim(-max(-abs(min(MinPrincStrain))*1.05, 0.37), max(abs(min(MinPrincStrain))*1.05, 0.37))
    plt.ylim(0, max(abs(max(MaxPrincStrain))*1.05, 0.39))
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))

    ax.annotate('Uniaxial \ncompression',xy=(0.145, 0.475), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Pure \nshear',xy=(0.13, 0.86), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Uniaxial \ntension',xy=(0.31, 0.868), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Plane strain \ntension',xy=(0.537, 0.87), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Equibiaxial \ntension',xy=(0.86, 0.7), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('FLC',xy=(0.68, 0.65), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))        
    ax.annotate('a',xy=(0.02, 0.96), xycoords='figure fraction', fontsize=28, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))    

    fig2.savefig(cwd + '\\Results\\FEA_stress_strain_plots\\' + test_config + '\\' + 'DIC_PrincStrain_' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)


    if show_plots == 'yes':
        plt.show()
    else:
        pass
        
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()