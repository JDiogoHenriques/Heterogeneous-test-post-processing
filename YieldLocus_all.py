"""
    description
   Plots principal stress and strain components

    INPUT: .csv file of the test with the LE, SDV1, S and EVOL field output request 
    OUTPUT: Plots the principal stress and strain components


    Run by writting on the command line: "python AbaqusPrincipalStressStrain_plots.py"
"""

# IMPORT PACKAGES -------------------------------------------------------------
from cmath import pi
from PythonScript_main import *
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from collections import Counter
import math
from scipy.interpolate import CubicSpline

#  --------------------------------------------------------------- method - MAIN
def main():

    cwd = os.getcwd()
    path_to_save_Results = os.path.join(cwd,'Results\\FEA_results\\')
    res_file_name = 'Load_'

    #Yield Locus processing data
    F = Ref_param[0]
    H = Ref_param[1]
    G = 1-H
    N = Ref_param[2]
    K = Ref_param[3]
    eps0 = Ref_param[4]
    n_swift = Ref_param[5]

    #Reading the number of data points per load step
    with open(path_to_save_Results + 'Strain_Stress_'+ test_config + '.csv', 'r') as csv_file:
        reader_aux = list(csv.reader(csv_file, delimiter=','))
        n_points_load = len(reader_aux[2:-1])+1

    #Importing and data management


    files = [filename for filename in os.listdir(path_to_save_Results) if
             filename.startswith(res_file_name) and filename.endswith("_AllSteps.csv")]
    if not files:
        print('ERROR: No .csv files with FEA results.')
        exit()
    
    all_results = []
    i = 0
    for file_name in files:  # read each file refering to a test
        with open(path_to_save_Results + file_name, 'r', errors='ignore') as infile:
            reader = list(csv.reader(infile, delimiter=','))
            all_results.append(np.array(reader[2:]))
            if i==0:
                n_points=np.array([len(all_results[i])])  # number of data points of first test
            else:
                n_points = np.append(n_points,len(all_results[i]))  # number of data points for each test
            i = i + 1  # counter

    EqPlastStrain = np.zeros((i,n_points[0],1))
    StressXX = np.zeros((i,n_points[0],1))
    StressYY  = np.zeros((i,n_points[0],1))
    StressXY = np.zeros((i,n_points[0],1))
    SigmaY = np.zeros((i,n_points[0],1))
    ElasticPoints = []
    ElastStressXX_norm = []
    ElastStressYY_norm = []
    PlastStressXX_norm = []
    PlastStressYY_norm = []

    for test in range(0, i):
        for line in range(0,n_points[0]):
            EqPlastStrain[test,line,0] = all_results[test][line][6]
            StressXX[test,line,0] = all_results[test][line][3]
            StressYY[test,line,0] = all_results[test][line][4]
            StressXY[test,line,0] = all_results[test][line][5]
            SigmaY[test,line,0] = K*((eps0+EqPlastStrain[test][line])**n_swift)
        ElasticPoints.append(np.array([x for x in range(len(EqPlastStrain[test])) if EqPlastStrain[test][x]==0]))
          
    StressXX_norm = StressXX/SigmaY
    StressYY_norm = StressYY/SigmaY
    StressXY_norm = StressXY/SigmaY

    x_min = -10
    x_max = 10
    y_min = -10
    y_max = 10
    x1, x2 = np.meshgrid(np.arange(x_min,x_max, 0.1), np.arange(y_min,y_max, 0.1))

    Hill_0 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0**2)
    Hill_02 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0.2**2)
    Hill_04 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0.4**2)
    Hill_06 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0.6**2) 

    for test_ind in range(0, i):
        if ElasticPoints[test_ind].size>0:
            ElastStressXX_norm.append(np.array(StressXX_norm[test_ind][ElasticPoints[test_ind]]))
            ElastStressYY_norm.append(np.array(StressYY_norm[test_ind][ElasticPoints[test_ind]]))

            PlastStressXX_norm.append(np.delete(StressXX_norm[test_ind], ElasticPoints[test_ind], 0))
            PlastStressYY_norm.append(np.delete(StressYY_norm[test_ind], ElasticPoints[test_ind], 0))
        else:
            pass 
# ------------- graph properties ----------------
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

    box = dict(boxstyle='Square,pad=0.05',fc='w',ec='w')
    al = 'center'
    al2 = 'left'
    fnt = 12
# ------------- Yield Locus - Material direction ------------- 

    colors = ['red','green','blue']

# ----------------------Fig 1 ------------------------------
    fig1 = plt.figure(figsize=(14,9))
    ax = plt.axes([0.22, 0.24, .86, .83])
    for j in range(len(available_test_configs)):

        if ElasticPoints[j].size>0:
            plastic = plt.scatter(PlastStressXX_norm[j],PlastStressYY_norm[j], c = colors[j], s=15,label= available_test_configs[j] + r'$\degree$ RD', zorder=2)
        else:
            plt.scatter(StressXX_norm[j],StressYY_norm[j], c = colors[j], s=15, zorder=2, label=available_test_configs[j] + r'$\degree$ RD')

    ElastStressXX_norm_aux = np.vstack([ElastStressXX_norm[i] for i in range(len(available_test_configs))])
    ElastStressYY_norm_aux = np.vstack([ElastStressYY_norm[i] for i in range(len(available_test_configs))])
    elastic = plt.scatter(ElastStressXX_norm_aux,ElastStressYY_norm_aux, c = 'grey', s=15,label=r'Elastic',zorder=1)

    YF0 = plt.contour(x1,x2,Hill_0,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF02 = plt.contour(x1,x2,Hill_02,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF04 = plt.contour(x1,x2,Hill_04,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF06 = plt.contour(x1,x2,Hill_06,[1], colors='k',linestyles='dashdot', linewidths=1.2)   

    plt.plot([0,0],[-1.5,1.5],'k--',lw=0.7)
    plt.plot([-1.5,1.5],[0,0],'k--',lw=0.7)

    #plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
    plt.xlabel('$\sigma_{11}/\sigma_y$')
    plt.ylabel('$\sigma_{22}/\sigma_y$')
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))

    name = '$\sigma_{12}/\sigma_y=0$, \nincrements of 0.2'
    plt.text(1.25,1.16,name,ha=al,va=al,ma=al2,bbox=box,fontsize=fnt)

    leg = ax.legend(loc='lower right',frameon=False, ncol=1, fontsize=fnt, markerscale=3.0,handletextpad=0.1,bbox_to_anchor=(0.1,0.01),
            handleheight=1.5, labelspacing=0.5)


    fig1.savefig(cwd + '\\Results\\YieldLocus\\' + test_config + '_YieldSurface_MaterialDir_Load' + available_test_configs[0] + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig1)

    # ----------------------Fig 2 ------------------------------

    fig2 = plt.figure(figsize=(14,9))
    ax = plt.axes([0.22, 0.24, .86, .83])
    for j in range(len(available_test_configs)):

        if ElasticPoints[j].size>0:
            plastic = plt.scatter(PlastStressXX_norm[len(available_test_configs)+j],PlastStressYY_norm[len(available_test_configs)+j], c = colors[j], s=15,label= available_test_configs[j] + r'$\degree$ RD', zorder=2)
        else:
            plt.scatter(StressXX_norm[len(available_test_configs)+j],StressYY_norm[len(available_test_configs)+j], c = colors[j], s=15, zorder=2, label=available_test_configs[j] + r'$\degree$ RD')

    ElastStressXX_norm_aux = np.vstack([ElastStressXX_norm[len(available_test_configs)+i] for i in range(len(available_test_configs))])
    ElastStressYY_norm_aux = np.vstack([ElastStressYY_norm[len(available_test_configs)+i] for i in range(len(available_test_configs))])
    elastic = plt.scatter(ElastStressXX_norm_aux,ElastStressYY_norm_aux, c = 'grey', s=15,label=r'Elastic',zorder=1)

    YF0 = plt.contour(x1,x2,Hill_0,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF02 = plt.contour(x1,x2,Hill_02,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF04 = plt.contour(x1,x2,Hill_04,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF06 = plt.contour(x1,x2,Hill_06,[1], colors='k',linestyles='dashdot', linewidths=1.2)   

    plt.plot([0,0],[-1.5,1.5],'k--',lw=0.7)
    plt.plot([-1.5,1.5],[0,0],'k--',lw=0.7)

    #plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
    plt.xlabel('$\sigma_{11}/\sigma_y$')
    plt.ylabel('$\sigma_{22}/\sigma_y$')
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))

    name = '$\sigma_{12}/\sigma_y=0$, \nincrements of 0.2'
    plt.text(1.25,1.16,name,ha=al,va=al,ma=al2,bbox=box,fontsize=fnt)

    leg = ax.legend(loc='lower right',frameon=False, ncol=1, fontsize=fnt, markerscale=3.0,handletextpad=0.1,bbox_to_anchor=(0.1,0.01),
            handleheight=1.5, labelspacing=0.5)


    fig2.savefig(cwd + '\\Results\\YieldLocus\\' + test_config + '_YieldSurface_MaterialDir_Load' + available_test_configs[1] + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig2)

    # ----------------------Fig 3 ------------------------------

    fig3 = plt.figure(figsize=(14,9))
    ax = plt.axes([0.22, 0.24, .86, .83])
    for j in range(len(available_test_configs)):

        if ElasticPoints[j].size>0:
            plastic = plt.scatter(PlastStressXX_norm[2*len(available_test_configs)+j],PlastStressYY_norm[2*len(available_test_configs)+j], c = colors[j], s=15,label= available_test_configs[j] + r'$\degree$ RD', zorder=2)
        else:
            plt.scatter(StressXX_norm[2*len(available_test_configs)+j],StressYY_norm[2*len(available_test_configs)+j], c = colors[j], s=15, zorder=2, label=available_test_configs[j] + r'$\degree$ RD')

    ElastStressXX_norm_aux = np.vstack([ElastStressXX_norm[2*len(available_test_configs)+i] for i in range(len(available_test_configs))])
    ElastStressYY_norm_aux = np.vstack([ElastStressYY_norm[2*len(available_test_configs)+i] for i in range(len(available_test_configs))])
    elastic = plt.scatter(ElastStressXX_norm_aux,ElastStressYY_norm_aux, c = 'grey', s=15,label=r'Elastic',zorder=1)

    YF0 = plt.contour(x1,x2,Hill_0,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF02 = plt.contour(x1,x2,Hill_02,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF04 = plt.contour(x1,x2,Hill_04,[1], colors='k',linestyles='dashdot', linewidths=1.2)
    YF06 = plt.contour(x1,x2,Hill_06,[1], colors='k',linestyles='dashdot', linewidths=1.2)   

    plt.plot([0,0],[-1.5,1.5],'k--',lw=0.7)
    plt.plot([-1.5,1.5],[0,0],'k--',lw=0.7)

    #plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
    plt.xlabel('$\sigma_{11}/\sigma_y$')
    plt.ylabel('$\sigma_{22}/\sigma_y$')
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))

    name = '$\sigma_{12}/\sigma_y=0$, \nincrements of 0.2'
    plt.text(1.25,1.16,name,ha=al,va=al,ma=al2,bbox=box,fontsize=fnt)

    leg = ax.legend(loc='lower right',frameon=False, ncol=1, fontsize=fnt, markerscale=3.0,handletextpad=0.1,bbox_to_anchor=(0.1,0.01),
            handleheight=1.5, labelspacing=0.5)


    fig3.savefig(cwd + '\\Results\\YieldLocus\\' + test_config + '_YieldSurface_MaterialDir_Load' + available_test_configs[2] + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig3)


    if show_plots == 'yes':
        plt.show()
    else:
        pass
        
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()