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
from matplotlib.ticker import FormatStrFormatter


#  --------------------------------------------------------------- method - MAIN
def main():

    cwd = os.getcwd()
    cwd_test = os.path.join(cwd,test_config,)

    num_ID_runs = 4
    num_iterations = 300
    num_stages = 40
    #path_DIC = os.path.join(cwd_test,'MatchID\\DIC\\')

    #Importing and data management
    Residual = np.zeros((num_iterations,num_ID_runs))
    Int_VirtualWork = np.zeros((num_stages,num_ID_runs))
    Ext_VirtualWork = np.zeros((num_stages,num_ID_runs))

    for i in range(num_ID_runs):
        path_to_Results = os.path.join(cwd_test,'MatchID\\VFM\\Results\\IdRun'+ str(i+1) +'\\')

        Residual_aux = []
        with open(path_to_Results + 'Residual.csv', 'r') as csv_file:
            reader = list(csv.reader(csv_file, delimiter=';'))
            n_points = len(reader)
            Residual_aux = np.array(reader[2:n_points])

        Residual_aux2 = Residual_aux[:,1]
        Residual_aux2 = np.array([float(x) for x in Residual_aux2])
    
        Residual[:,i] = Residual_aux2

        VirtualWork_aux = []
        with open(path_to_Results + 'VirtualWork.csv', 'r') as csv_file:
            reader = list(csv.reader(csv_file, delimiter=';'))
            n_points = len(reader)
            VirtualWork_aux = np.array(reader[2:n_points])

        Ext_VirtualWork_aux2 = np.array([x for x in VirtualWork_aux[0:num_stages]])
        Ext_VirtualWork_aux3 = Ext_VirtualWork_aux2[:,1]
        Ext_VirtualWork_aux3 = np.array([float(x) for x in Ext_VirtualWork_aux3])
        Ext_VirtualWork[:,i] = Ext_VirtualWork_aux3

        Int_VirtualWork_aux2 = np.array([x for x in VirtualWork_aux[num_stages+2:-2]])
        Int_VirtualWork_aux3 = np.array([x for x in Int_VirtualWork_aux2[0:num_stages]])
        Int_VirtualWork_aux3 = Int_VirtualWork_aux3[:,1]
        Int_VirtualWork_aux3 = np.array([float(x) for x in Int_VirtualWork_aux3])
        Int_VirtualWork[:,i] = Int_VirtualWork_aux3

    cwd = os.getcwd()
    path_to_save_Results = os.path.join(cwd,'Results\\FEA_results\\')

    #Reading the number of data points per load step
    with open(path_to_save_Results + 'Strain_Stress_'+ test_config + '.csv', 'r') as csv_file:
        reader_aux = list(csv.reader(csv_file, delimiter=','))
        n_points_load = len(reader_aux[2:-1])+1

    #Importing and data management
    all_results = []
    with open(path_to_save_Results + test_config + '_AllSteps.csv', 'r') as csv_file:
        reader = list(csv.reader(csv_file, delimiter=','))
        n_points = len(reader)
        all_results = np.array(reader[2:n_points])

    EqPlastStrain = all_results[n_points_load:-1,6]
    StressXX = all_results[n_points_load:-1,3]
    StressYY = all_results[n_points_load:-1,4]
    StressXY = all_results[n_points_load:-1,5]
    MinPrincStress = all_results[n_points_load:-1,7]
    MaxPrincStress = all_results[n_points_load:-1,8]

    EqPlastStrain = np.array([float(x) for x in EqPlastStrain])
    StressXX = np.array([float(x) for x in StressXX])
    StressYY = np.array([float(x) for x in StressYY])
    StressXY = np.array([float(x) for x in StressXY])
    MinPrincStress = np.array([float(x) for x in MinPrincStress])
    MaxPrincStress = np.array([float(x) for x in MaxPrincStress])

    #Yield Locus processing data
    F = Ref_param[0]
    H = Ref_param[1]
    G = 1-H
    N = Ref_param[2]
    K = Ref_param[3]
    eps0 = Ref_param[4]
    n_swift = Ref_param[5]

    SigmaY = K*((eps0+EqPlastStrain)**n_swift)

    StressXX_norm = StressXX/SigmaY
    StressYY_norm = StressYY/SigmaY
    StressXY_norm = StressXY/SigmaY
    MinPrincStress_norm = MinPrincStress/SigmaY
    MaxPrincStress_norm = MaxPrincStress/SigmaY

    x_min = -10
    x_max = 10
    y_min = -10
    y_max = 10
    x1, x2 = np.meshgrid(np.arange(x_min,x_max, 0.1), np.arange(y_min,y_max, 0.1))

    Hill_0 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0**2)
    Hill_02 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0.25**2)
    Hill_04 = H*(x1-x2)**2+G*(x1**2)+F*(x2**2)+2*N*(0.5**2)

    #Id Run 1
    F_1 = 0.4185
    H_1 = 0.3303
    G_1 = 1-H
    N_1 = 1.2970
    K_1 = 1054.00
    eps0_1 = 0.00534
    n_swift_1 = 0.199

    Hill_0_1 = H_1*(x1-x2)**2+G_1*(x1**2)+F_1*(x2**2)+2*N_1*(0**2)
    Hill_02_1 = H_1*(x1-x2)**2+G_1*(x1**2)+F_1*(x2**2)+2*N_1*(0.25**2)
    Hill_04_1 = H_1*(x1-x2)**2+G_1*(x1**2)+F_1*(x2**2)+2*N_1*(0.5**2)

    #Id Run 2
    F_2 = 0.3925
    H_2 = 0.3301
    G_2 = 1-H
    N_2 = 1.25
    K_2 = 1036
    eps0_2 = 0.00559
    n_swift_2 = 0.201

    Hill_0_2 = H_2*(x1-x2)**2+G_2*(x1**2)+F_2*(x2**2)+2*N_2*(0**2)
    Hill_02_2 = H_2*(x1-x2)**2+G_2*(x1**2)+F_2*(x2**2)+2*N_2*(0.25**2)
    Hill_04_2 = H_2*(x1-x2)**2+G_2*(x1**2)+F_2*(x2**2)+2*N_2*(0.5**2)

    #Id Run 3
    F_3 = 0.4872
    H_3 = 0.3296
    G_3 = 1-H
    N_3 = 1.4280
    K_3 = 1086.00
    eps0_3 = 0.00371
    n_swift_3 = 0.1875

    Hill_0_3 = H_3*(x1-x2)**2+G_3*(x1**2)+F_3*(x2**2)+2*N_3*(0**2)
    Hill_02_3 = H_3*(x1-x2)**2+G_3*(x1**2)+F_3*(x2**2)+2*N_3*(0.25**2)
    Hill_04_3 = H_3*(x1-x2)**2+G_3*(x1**2)+F_3*(x2**2)+2*N_3*(0.5**2)
    SigmaY_3 = K_3*((eps0_3+EqPlastStrain)**n_swift_3)

    #Id Run 4
    F_4 = 0.4872
    H_4 = 0.3296
    G_4 = 1-H
    N_4 = 1.2430
    K_4 = 1043.00
    eps0_4 = 0.00700
    n_swift_4 = 0.209

    Hill_0_4 = H_4*(x1-x2)**2+G_4*(x1**2)+F_4*(x2**2)+2*N_4*(0**2)
    Hill_02_4 = H_4*(x1-x2)**2+G_4*(x1**2)+F_4*(x2**2)+2*N_4*(0.25**2)
    Hill_04_4 = H_4*(x1-x2)**2+G_4*(x1**2)+F_4*(x2**2)+2*N_4*(0.5**2)

    #Hardening curve processing
    EqPlasticStrain_hardening = np.unique(sorted(EqPlastStrain))
    SigmaY_hardening = K*((eps0+EqPlasticStrain_hardening)**n_swift)
    SigmaY_hardening_1 = K_1*((eps0_1+EqPlasticStrain_hardening)**n_swift_1)
    SigmaY_hardening_2 = K_2*((eps0_2+EqPlasticStrain_hardening)**n_swift_2)
    SigmaY_hardening_3 = K_3*((eps0_3+EqPlasticStrain_hardening)**n_swift_3)
    SigmaY_hardening_4 = K_4*((eps0_4+EqPlasticStrain_hardening)**n_swift_4)

    #Checking if there are points in elasticity
    ElasticPoints = np.array([x for x in range(len(EqPlastStrain)) if EqPlastStrain[x]==0])

    if ElasticPoints.size>0:
        ElastStressXX_norm = np.array(StressXX_norm[ElasticPoints])
        ElastStressYY_norm = np.array(StressYY_norm[ElasticPoints])
        ElastStressXY_norm = np.array(StressXY_norm[ElasticPoints])

        PlastStressXX_norm = np.delete(StressXX_norm, ElasticPoints, 0)
        PlastStressYY_norm = np.delete(StressYY_norm, ElasticPoints, 0)
        PlastStressXY_norm = np.delete(StressXY_norm, ElasticPoints, 0)

    else:
        pass    

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
    box = dict(boxstyle='Square,pad=0.05',fc='w',ec='w')
    al = 'center'
    al2 = 'left'
    fnt = 18

# -------------Residual graph------------- 
    fig1 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(Residual[:,0], color='red', label="Id. Run 1", linewidth=2.5) # Id. Run 1
    plt.plot(Residual[:,1], color='black', label="Id. Run 2", linewidth=2.5) # Id. Run 2
    plt.plot(Residual[:,2], color='blue', label="Id. Run 3", linewidth=2.5) # Id. Run 3
    plt.plot(Residual[:,3], color='grey', label="Id. Run 4", linewidth=2.5) # Id. Run 4

    leg = ax.legend(loc='upper right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Principal strain standard deviation")
    plt.xlabel("Iterations", fontsize=FS, labelpad=10)
    plt.ylabel("Residual [$\mathrm{N}^{2}.\mathrm{mm}^{2}$]", fontsize=FS, labelpad=10)
       
    ax.set_yscale('symlog')
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    plt.xlim(0, num_iterations+1)
    plt.ylim(-0.001, max(Residual[0]))
    #ax.yaxis.set_major_formatter(mticker.ScalarFormatter())

    fig1.savefig(cwd + '\\Results\\Identification\\' + 'Residual.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)

# -------------Virtual Work graph------------- 
    fig2 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.12, 0.14, .86, .83])

    plt.plot(range(0,num_stages+1),np.append(0,Ext_VirtualWork[:,0]), color='black', label="Ext. virtual work", linewidth=2.5) # Ext. virtual work

    plt.plot(range(0,num_stages+1),np.append(0,Int_VirtualWork[:,0]), color='red', ls='--', marker='v', zorder=10, markersize=8, clip_on=False, label="Int. virtual work - Id. Run 1", linewidth=2.5) # Id. Run 1
    plt.plot(range(0,num_stages+1),np.append(0,Int_VirtualWork[:,0]), color='black', ls='--', marker='o', zorder=10, markersize=8, clip_on=False, label="Int. virtual work - Id. Run 2", linewidth=2.5) # Id. Run 2
    plt.plot(range(0,num_stages+1),np.append(0,Int_VirtualWork[:,0]), color='blue', ls='--', marker='^', zorder=10, markersize=8, clip_on=False, label="Int. virtual work - Id. Run 3", linewidth=2.5) # Id. Run 3
    plt.plot(range(0,num_stages+1),np.append(0,Int_VirtualWork[:,0]), color='grey', ls='--', marker='.', zorder=10, markersize=8, clip_on=False, label="Int. virtual work - Id. Run 4", linewidth=2.5) # Id. Run 4

    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Principal strain standard deviation")
    plt.xlabel("Load step", fontsize=FS, labelpad=10)
    plt.ylabel("Virtual work [N.mm]", fontsize=FS, labelpad=10)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    plt.xlim(0, num_stages)
    plt.ylim(0, max(Int_VirtualWork[-1])*1.05)
    plt.locator_params(axis='x', nbins=4)
    plt.locator_params(axis='y', nbins=5)
    fig2.savefig(cwd + '\\Results\\Identification\\' + 'VirtualWork.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)


# ------------- Yield Locus - Material direction ------------- 

    fig3 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.22, 0.24, .86, .83])

    colors = EqPlastStrain
    if ElasticPoints.size>0:
        colors = np.delete(colors, ElasticPoints, 0)
        elastic = plt.scatter(ElastStressXX_norm,ElastStressYY_norm, c = 'grey', s=15,label=r'Elastic')
        plastic = plt.scatter(PlastStressXX_norm,PlastStressYY_norm, c = colors, cmap='jet', s=15,label=r'Plastic', alpha= 1)
    else:
        plt.scatter(StressXX_norm,StressYY_norm, c = colors, cmap='jet', s=15,label=r'Plastic')

    YF0 = plt.contour(x1,x2,Hill_0,[1], colors='k',linestyles='dashdot', linewidths=1.5)
    YF02 = plt.contour(x1,x2,Hill_02,[1], colors='k',linestyles='dashdot', linewidths=1.5)
    YF04 = plt.contour(x1,x2,Hill_04,[1], colors='k',linestyles='dashdot', linewidths=1.5)

    YF0_1 = plt.contour(x1,x2,Hill_0_1,[1], colors='red', linewidths=1.2)
    YF02_1 = plt.contour(x1,x2,Hill_02_1,[1], colors='red', linewidths=1.2)
    YF04_1 = plt.contour(x1,x2,Hill_04_1,[1], colors='red', linewidths=1.2)

    YF0_2 = plt.contour(x1,x2,Hill_0_2,[1], colors='g', linewidths=1.2)
    YF02_2 = plt.contour(x1,x2,Hill_02_2,[1], colors='g', linewidths=1.2)
    YF04_2 = plt.contour(x1,x2,Hill_04_2,[1], colors='g', linewidths=1.2)

    YF0_3 = plt.contour(x1,x2,Hill_0_3,[1], colors='blue', linewidths=1.5)
    YF02_3 = plt.contour(x1,x2,Hill_02_3,[1], colors='blue', linewidths=1.5)
    YF04_3 = plt.contour(x1,x2,Hill_04_3,[1], colors='blue', linewidths=1.5)

    YF0_4 = plt.contour(x1,x2,Hill_0_4,[1], colors='grey', linewidths=1.2)
    YF02_4 = plt.contour(x1,x2,Hill_02_4,[1], colors='grey', linewidths=1.2)
    YF04_4 = plt.contour(x1,x2,Hill_04_4,[1], colors='grey', linewidths=1.2)

    plt.plot([0,0],[-1.5,1.5],'k--',lw=0.7)
    plt.plot([-1.5,1.5],[0,0],'k--',lw=0.7)


    plt.xlabel('$\sigma_{11}/\sigma_y$')
    plt.ylabel('$\sigma_{22}/\sigma_y$')
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))

    name = '$\sigma_{12}/\sigma_y=0$, \nincrements of 0.25'
    plt.text(1.09,1.37,name,ha=al,va=al,ma=al2,bbox=box,fontsize=fnt-3.3)

    h0,l0 = YF0.legend_elements()
    h1,l1 = YF0_1.legend_elements()
    h2,l2 = YF0_2.legend_elements()
    h3,l3 = YF0_3.legend_elements()
    h4,l4 = YF0_4.legend_elements()
    
    leg = plt.legend([elastic,plastic,h0[0],h1[0],h2[0],h3[0],h4[0]],['Elastic','Plastic','Reference','Id. Run 1','Id. Run 2','Id. Run 3','Id. Run 4'],loc='lower right',frameon=False, ncol=1, fontsize=fnt-3.5, markerscale=3.0,handletextpad=0.1,
                handleheight=1.5, labelspacing=0.5)

    for l in leg.get_lines():
        l.set_alpha(1)

    cjet = plt.cm.get_cmap('jet', 12)
    sm = plt.cm.ScalarMappable(cmap=cjet,norm=plt.Normalize(vmin=0, vmax=max(EqPlastStrain)))
    sm._A = []
    tks = np.linspace(0.0,max(EqPlastStrain),9)
    clb = plt.colorbar(sm,ticks=tks, ax=ax)
    clb.solids.set_rasterized(True)
    clb.set_label(r'$\bar \epsilon ^\mathrm{P}$',labelpad=-65,y=1.07, rotation=0)
    clb.set_ticklabels(['%.2f'%round(i,2) for i in tks])
    plt.title(" ")
    clb.solids.set(alpha=1)

    plt.clim(0,max(abs(EqPlastStrain)))
    fig3.savefig(cwd + '\\Results\\Identification\\' + 'YieldLocus.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig1)

# ------------- Swift Hardening law ------------- 

    fig4 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.22, 0.24, .86, .83])

    plt.plot(EqPlasticStrain_hardening,SigmaY_hardening, color='black',ls='--', label="Reference", linewidth=2.5) 

    plt.plot(EqPlasticStrain_hardening,SigmaY_hardening_1, color='red', label="Id. Run 1", linewidth=2.5) # Id. Run 1
    plt.plot(EqPlasticStrain_hardening,SigmaY_hardening_2, color='black', label="Id. Run 2", linewidth=2.5) # Id. Run 2
    plt.plot(EqPlasticStrain_hardening,SigmaY_hardening_3, color='blue', label="Id. Run 3", linewidth=2.5) # Id. Run 3
    plt.plot(EqPlasticStrain_hardening,SigmaY_hardening_4, color='grey', label="Id. Run 4", linewidth=2.5) # Id. Run 4

    leg = ax.legend(loc='lower right', fontsize=FS-2, frameon=True, ncol=1, handleheight=1.5, labelspacing=0.5)
    #plt.title("Principal strain standard deviation")
    plt.xlabel(r'$\bar \epsilon ^\mathrm{P}$', fontsize=FS, labelpad=10)
    plt.ylabel(r'$\sigma_y$ [MPa]', fontsize=FS, labelpad=10)

    ax.tick_params(axis='x', colors=(0,0,0))
    ax.tick_params(axis='y', colors=(0,0,0))
    plt.xlim(0, max(EqPlasticStrain_hardening)*1.01)
    plt.ylim(0, max(SigmaY_hardening)*1.2)
    plt.locator_params(axis='x', nbins=6)
    plt.locator_params(axis='y', nbins=6)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    fig4.savefig(cwd + '\\Results\\Identification\\' + 'Stress_Strain_HardeningPlot.jpg', dpi=600, bbox_inches='tight', pad_inches=.05)


    if show_plots == 'yes':
        plt.show()
    else:
        pass
        
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()