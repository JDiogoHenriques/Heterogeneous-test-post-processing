"""
"""

# IMPORT PACKAGES -------------------------------------------------------------
from PythonScript_main import *
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

#  --------------------------------------------------------------- method - MAIN
def main():

    Step = 40
    S_size= 23
    Transparency = 0.9
    Mark = 's'
    cwd = os.getcwd()
    path_to_save_Results = os.path.join(cwd,'Results\\Sensitivity_analysis\\')
    path_to_Results_FEA = os.path.join(cwd,'Results\\FEA_results\\')

    #Reading the number of data points
    with open(path_to_Results_FEA + 'Strain_Stress_'+ test_config + '.csv', 'r') as csv_file:
        reader = list(csv.reader(csv_file, delimiter=','))
        n_points = len(reader)-2

    #Reading the number of increments
    with open(path_to_save_Results + test_config + '_' + Param[0] + '_Perturbation.csv', 'r') as csv_file:
        reader = list(csv.reader(csv_file, delimiter=','))
        n_points_inc = len(reader)
        n_inc = int(n_points_inc/n_points)
    
    #Importing and data management - Perturbation
    Perturbation_all_results = np.zeros((n_points*n_inc, 9, len(Param))) 
    for i in range(len(Param)):
        with open(path_to_save_Results + test_config + '_' + Param[i] + '_Perturbation.csv', 'r') as csv_file:
            reader = list(csv.reader(csv_file, delimiter=','))
            n_points_pert = len(reader)
            Perturbation_all_results[:,:,i] = np.array(reader[2:n_points_pert])

    #Importing and data management - Reference
    Reference_all_results = np.zeros((n_points*n_inc, 6))
    with open(cwd + '\\Results\\FEA_results\\' + test_config + '_AllSteps.csv', 'r') as csv_file:
        reader = list(csv.reader(csv_file, delimiter=','))
        n_points_ref = len(reader)
        Reference_all_results = np.array(reader[2:n_points_ref], dtype=float)#.astype(np.float)

    #Calculating differences
    Strain_xx_diff = np.zeros((n_points*n_inc, len(Param)))
    Strain_yy_diff = np.zeros((n_points*n_inc, len(Param)))
    Strain_xy_diff = np.zeros((n_points*n_inc, len(Param)))
    Sigma_xx_diff = np.zeros((n_points*n_inc, len(Param)))
    Sigma_yy_diff = np.zeros((n_points*n_inc, len(Param)))
    Sigma_xy_diff = np.zeros((n_points*n_inc, len(Param)))

    for i in range(len(Param)):
        Strain_xx_diff[:,i] = Reference_all_results[:,0] - Perturbation_all_results[:,0,i]
        Strain_yy_diff[:,i] = Reference_all_results[:,1] - Perturbation_all_results[:,1,i]
        Strain_xy_diff[:,i] = Reference_all_results[:,2] - Perturbation_all_results[:,2,i]

        Sigma_xx_diff[:,i] = Reference_all_results[:,3] - Perturbation_all_results[:,3,i]
        Sigma_yy_diff[:,i] = Reference_all_results[:,4] - Perturbation_all_results[:,4,i]
        Sigma_xy_diff[:,i] = Reference_all_results[:,5] - Perturbation_all_results[:,5,i]

    Strain_xx_diff_byinc = np.zeros((n_inc, len(Param)))
    Strain_yy_diff_byinc = np.zeros((n_inc, len(Param)))
    Strain_xy_diff_byinc = np.zeros((n_inc, len(Param)))

    for k in range(n_inc):
        for j in range(len(Param)):
            Strain_xx_diff_byinc [k,j] = np.sqrt(np.sum((Strain_xx_diff[((k)*n_points):((k+1)*(n_points))][:,j]/np.amax(abs(Reference_all_results),axis=0)[0])**2))
            Strain_yy_diff_byinc [k,j] = np.sqrt(np.sum((Strain_yy_diff[((k)*n_points):((k+1)*(n_points))][:,j]/np.amax(abs(Reference_all_results),axis=0)[1])**2))
            Strain_xy_diff_byinc [k,j] = np.sqrt(np.sum((Strain_xy_diff[((k)*n_points):((k+1)*(n_points))][:,j]/np.amax(abs(Reference_all_results),axis=0)[2])**2))
    
    F_SensRank = (1/(3*n_inc*n_points))*(Strain_xx_diff_byinc[:,0].sum()+Strain_yy_diff_byinc[:,0].sum()+Strain_xy_diff_byinc[:,0].sum())
    H_SensRank = (1/(3*n_inc*n_points))*(Strain_xx_diff_byinc[:,1].sum()+Strain_yy_diff_byinc[:,1].sum()+Strain_xy_diff_byinc[:,1].sum())
    N_SensRank = (1/(3*n_inc*n_points))*(Strain_xx_diff_byinc[:,2].sum()+Strain_yy_diff_byinc[:,2].sum()+Strain_xy_diff_byinc[:,2].sum())
    K_SensRank = (1/(3*n_inc*n_points))*(Strain_xx_diff_byinc[:,3].sum()+Strain_yy_diff_byinc[:,3].sum()+Strain_xy_diff_byinc[:,3].sum())
    eps_0_SensRank = (1/(3*n_inc*n_points))*(Strain_xx_diff_byinc[:,4].sum()+Strain_yy_diff_byinc[:,4].sum()+Strain_xy_diff_byinc[:,4].sum())    
    n_swift_SensRank = (1/(3*n_inc*n_points))*(Strain_xx_diff_byinc[:,5].sum()+Strain_yy_diff_byinc[:,5].sum()+Strain_xy_diff_byinc[:,5].sum())
  
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

# -------------Sensitivity bar plot------------- 
    fig1 = plt.figure(figsize=(11,9))
    ax = plt.axes([0.12, 0.14, .86, .83])

    ax.barh(Param, [F_SensRank,H_SensRank,N_SensRank,K_SensRank,eps_0_SensRank,n_swift_SensRank])

    #leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1,
    #            handleheight=1.5, labelspacing=0.5)
    plt.title("Test configuration: " + test_config + " | Hill48 + Swift law")
    plt.xlabel("Parameters sensitivity", fontsize=FS, labelpad=10)
    plt.ylabel("Parameters", fontsize=FS, labelpad=10)
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    fig1.savefig(cwd + '\\Results\\Sensitivity_analysis\\Sensitivity_barPlot.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)
    plt.close(fig1)

    print("Plotting difference fields for step: " + str(Step))
# -------------Diff EpsX------------- 
    for i in range(len(Param)):
        fig2 = plt.figure(figsize=(11,7))
        ax = plt.axes([0.01, 0.01, .9, .95])

        plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Strain_xx_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'EpsXX_diff')
    
        plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.xlim(-20,80)
        plt.ylim(-20,80)

        cbar = plt.colorbar()

        plt.clim(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
        #sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
        tks = np.linspace(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\varepsilon_{xx}$',labelpad=-40,y=1.1, rotation=0)
        #cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
        plt.axis('off')
        cbar.solids.set(alpha=1)
        fig2.savefig(cwd + '\\Results\\Sensitivity_analysis\\EpsXXdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
        plt.close(fig2)
# -------------Diff Epsy------------- 
    for i in range(len(Param)):
        fig3 = plt.figure(figsize=(11,7))
        ax = plt.axes([0.01, 0.01, .9, .95])

        plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Strain_yy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'EpsYY_diff')
    
        plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.xlim(-20,80)
        plt.ylim(-20,80)

        cbar = plt.colorbar()

        plt.clim(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
        #sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
        tks = np.linspace(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)       
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\varepsilon_{yy}$',labelpad=-40,y=1.1, rotation=0)
        #cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
        plt.axis('off')
        cbar.solids.set(alpha=1)
        fig3.savefig(cwd + '\\Results\\Sensitivity_analysis\\EpsYYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
        plt.close(fig3)
# -------------Diff Epsxy------------- 
    for i in range(len(Param)):
        fig4 = plt.figure(figsize=(11,7))
        ax = plt.axes([0.01, 0.01, .9, .95])

        plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Strain_xy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'EpsYY_diff')
    
        plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.xlim(-20,80)
        plt.ylim(-20,80)

        cbar = plt.colorbar()

        plt.clim(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
        #sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
        tks = np.linspace(0,max(abs(Strain_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\varepsilon_{xy}$',labelpad=-40,y=1.1, rotation=0)
        #cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
        plt.axis('off')
        cbar.solids.set(alpha=1)
        fig4.savefig(cwd + '\\Results\\Sensitivity_analysis\\EpsXYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
        plt.close(fig4)
# -------------Diff SigmaX------------- 
    for i in range(len(Param)):
        fig5 = plt.figure(figsize=(11,7))
        ax = plt.axes([0.01, 0.01, .9, .95])

        plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Sigma_xx_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'SigmaXX_diff')
    
        plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.xlim(-20,80)
        plt.ylim(-20,80)

        cbar = plt.colorbar()

        plt.clim(0,max(abs(Sigma_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
        #sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
        tks = np.linspace(0,max(abs(Sigma_xx_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\sigma_{xx}$',labelpad=-40,y=1.1, rotation=0)
        #cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
        plt.axis('off')
        cbar.solids.set(alpha=1)
        fig5.savefig(cwd + '\\Results\\Sensitivity_analysis\\SigmaXXdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
        plt.close(fig5)
# -------------Diff SigmaY------------- 
    for i in range(len(Param)):
        fig6 = plt.figure(figsize=(11,7))
        ax = plt.axes([0.01, 0.01, .9, .95])

        plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Sigma_yy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'SigmaYY_diff')
    
        plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.xlim(-20,80)
        plt.ylim(-20,80)

        cbar = plt.colorbar()

        plt.clim(0,max(abs(Sigma_yy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
        #sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
        tks = np.linspace(0,max(abs(Sigma_yy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\sigma_{yy}$',labelpad=-40,y=1.1, rotation=0)
        #cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
        plt.axis('off')
        cbar.solids.set(alpha=1)
        fig6.savefig(cwd + '\\Results\\Sensitivity_analysis\\SigmaYYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
        plt.close(fig6)
# -------------Diff SigmaXY------------- 
    for i in range(len(Param)):
        fig7 = plt.figure(figsize=(11,7))
        ax = plt.axes([0.01, 0.01, .9, .95])

        plt.scatter(Perturbation_all_results[0:n_points,6,0],Perturbation_all_results[0:n_points,7,0], marker=Mark, c = abs(Sigma_xy_diff[n_points*(Step-1):Step*n_points,i]), cmap='jet', alpha =Transparency, s=S_size, label=r'SigmaXY_diff')
    
        plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.xlim(-20,80)
        plt.ylim(-20,80)

        cbar = plt.colorbar()

        plt.clim(0,max(abs(Sigma_xy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])))
        #sm1 = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=min(Strain_xx_diff.min(axis=1)[n_points*(n_inc-1):n_inc*n_points]),vmax=max(Strain_xx_diff.max(axis=1)[n_points*(n_inc-1):n_inc*n_points])))
        tks = np.linspace(0,max(abs(Sigma_xy_diff.max(axis=1)[n_points*(Step-1):Step*n_points])),13)
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\sigma_{xy}$',labelpad=-40,y=1.1, rotation=0)
        #cbar.set_ticklabels(['%.3f'%round(i,3) for i in tks])
        plt.axis('off')
        cbar.solids.set(alpha=1)
        fig7.savefig(cwd + '\\Results\\Sensitivity_analysis\\SigmaXYdiff_'+ Param[i] +'_Perturbation.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
        plt.close(fig7)

    if show_plots == 'yes':
        plt.show()
    else:
        pass
        
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()