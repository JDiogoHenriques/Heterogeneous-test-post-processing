"""
    description
   Plots principal stress and strain components

    INPUT: .csv file of the test with the LE, SDV1, S and EVOL field output request 
    OUTPUT: Plots the principal stress and strain components


    Run by writting on the command line: "python AbaqusPrincipalStressStrain_plots.py"
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

    #Importing and data management
    all_results = []
    with open(path_to_save_Results + 'Strain_Stress_'+ test_config + '.csv', 'r') as csv_file:
        reader = list(csv.reader(csv_file, delimiter=','))
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

#--------------- Stress processing ------------------------------
#Stress Uniaxial compression
    Stress_uni_comp_x = np.linspace(-10000,10000)
    Stress_uni_comp_y = Stress_uni_comp_x*0

#Stress pure shear
    Stress_pure_shear_x = np.linspace(-10000,10000)
    Stress_pure_shear_y = Stress_pure_shear_x/-1

#Stress uniaxial tension
    Stress_uni_tens_y = np.linspace(-10000,10000)
    Stress_uni_tens_x = Stress_uni_tens_y*0

#Stress equibiaxial tension
    Stress_equibiaxial_tens_x = np.linspace(-10000,10000)
    Stress_equibiaxial_tens_y = Stress_equibiaxial_tens_x

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
    fig1 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.22, 0.24, .86, .83])

    colors = EqPlastStrain
    if ElasticPoints.size>0:
        colors = np.delete(colors, ElasticPoints, 0)
        elastic = plt.scatter(ElastMinPrincStrain,ElastMaxPrincStrain, s=35, c = 'grey', marker = '+', label=r'Elastic')
        plastic = plt.scatter(PlastMinPrincStrain,PlastMaxPrincStrain, c = colors, cmap='jet', s=15, label=r'Plastic')
    else:
        elastic = plt.scatter([],[], s=35, c = 'grey', marker = '+', label=r'Elastic')
        plastic = plt.scatter(MinPrincStrain,MaxPrincStrain, c = colors, cmap='jet', s=15, label=r'Plastic')
            
    plt.plot(Strain_uni_comp_x, Strain_uni_comp_y, color='black', linewidth=0.8)
    plt.plot(Strain_pure_shear_x, Strain_pure_shear_y, color='black', linewidth=0.8)
    plt.plot(Strain_uni_tens_x,Strain_uni_tens_y, color='black', linewidth=0.8)
    plt.plot(Strain_planestra_tens_x,Strain_planestra_tens_y, color='black', linewidth=0.8)
    plt.plot(Strain_equibiaxial_tens_x,Strain_equibiaxial_tens_y, color='black', linewidth=0.8)
    plt.plot(FLC_Eps2,FLC_Eps1, color='grey', linewidth=0.7)

    leg = ax.legend(loc='lower right', fontsize=FS-6, frameon=True, ncol=1,
                handleheight=1.5, labelspacing=0.5)
    plt.title(" ")
    plt.xlabel(r'$\epsilon_2$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\epsilon_1$', fontsize=FS, labelpad=10)
    plt.xlim(-max(abs(min(MinPrincStrain))*1.05, 0.38), max(abs(min(MinPrincStrain))*1.05, 0.38))
    plt.ylim(-max(-abs(min(MaxPrincStrain))*1.05, 0), max(abs(max(MaxPrincStrain))*1.05, 0.49))
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    #plt.colorbar()

    ax.annotate('Uniaxial \ncompression',xy=(0.17, 0.42), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Pure \nshear',xy=(0.15, 0.72), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Uniaxial \ntension',xy=(0.29, 0.9), zorder=5, xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Plane strain \ntension',xy=(0.535, 0.87), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Equibiaxial \ntension',xy=(0.88, 0.63), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('FLC',xy=(0.73, 0.59), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))        
    #ax.annotate(r'$\bar{\epsilon}^\mathrm{p}$',xy=(0.855, 0.96), xycoords='figure fraction', fontsize=28, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))    
    #ax.annotate('a',xy=(0.08, 0.96), xycoords='figure fraction', fontsize=28, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))    

    fig1.savefig(cwd + '\\Results\\FEA_stress_strain_plots\\' + test_config + '\\' + 'PrincStrain_' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

# -------------Principal stress graph-------------     
    fig2 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    if ElasticPoints.size>0:
        elasticstress = plt.scatter(ElastMinPrincStress,ElastMaxPrincStress, s=35, c = 'grey', marker = '+', label=r'Elastic')
        plasticstrain = plt.scatter(PlastMinPrincStress,PlastMaxPrincStress, c = colors, cmap='jet', s=15, label=r'Plastic')
    else:
        elasticstress = plt.scatter([],[], s=35, c = 'grey', marker = '+', label=r'Elastic')
        plasticstrain = plt.scatter(MinPrincStress,MaxPrincStress, c = colors, cmap='jet', s=15, label=r'Plastic')

    plt.plot(Stress_uni_comp_x,Stress_uni_comp_y, color='black', linewidth=0.8)
    plt.plot(Stress_pure_shear_x,Stress_pure_shear_y, color='black', linewidth=0.8)
    plt.plot(Stress_uni_tens_x,Stress_uni_tens_y, color='black', linewidth=0.8)
    plt.plot(Stress_equibiaxial_tens_x,Stress_equibiaxial_tens_y, color='black', linewidth=0.8)

    leg = ax.legend(loc='lower right', fontsize=FS-6, frameon=True, ncol=1,
                handleheight=1.5, labelspacing=0.5)
    plt.title(" ")
    plt.xlabel(r'$\sigma_2$ [MPa]', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\sigma_1$ [MPa]', fontsize=FS, labelpad=10)
    plt.xlim(-max(abs(min(MinPrincStress))*1.1, 1200), max(abs(min(MinPrincStress))*1.1, 1200))
    plt.ylim(-max(-220+min(MinPrincStress), 300), max(abs(max(MaxPrincStress))*1.1, 1190))
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    #plt.colorbar()

    ax.annotate('Uniaxial \ncompression',xy=(0.195, 0.315), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Pure \nshear',xy=(0.22, 0.84), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Uniaxial \ntension',xy=(0.565, 0.91), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Equibiaxial \ntension',xy=(0.905, 0.795), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))
    ax.annotate('Equibiaxial \ncompression',xy=(0.485, 0.15), xycoords='figure fraction', fontsize=12, horizontalalignment='center', verticalalignment='center',bbox = dict(boxstyle ="round", fc ="1", ec="1"))    

    fig2.savefig(cwd + '\\Results\\FEA_stress_strain_plots\\' + test_config + '\\' + 'PrincStress_' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

    if show_plots == 'yes':
        plt.show()
    else:
        pass
        
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()