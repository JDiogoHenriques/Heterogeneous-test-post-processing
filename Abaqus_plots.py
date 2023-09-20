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

    #test_name = 'ArcanSpecimen'  # USE NAME OF TEST
    #test_config = 'Load_0_MatDir_0' #Load_90_MatDir_0 , Load_45_MatDir_0, Load_0_MatDir_0,...
    cwd = os.getcwd()
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
    IntPointCoordX = all_results[:,6]
    IntPointCoordY = all_results[:,7]
    StressXX = all_results[:,8]
    StressYY = all_results[:,9]
    StressXY = all_results[:,10]
    StressTriax = all_results[:,11]
    MisesStress = all_results[:,12]
    ThirdINVStress = all_results[:,13]

    MinPrincStrain = np.array([float(x) for x in MinPrincStrain])
    MaxPrincStrain = np.array([float(x) for x in MaxPrincStrain])
    EqPlastStrain = np.array([float(x) for x in EqPlastStrain])
    MinPrincStress = np.array([float(x) for x in MinPrincStress])
    MaxPrincStress = np.array([float(x) for x in MaxPrincStress])
    IntPointCoordX = np.array([float(x) for x in IntPointCoordX])
    IntPointCoordY = np.array([float(x) for x in IntPointCoordY])
    StressXX = np.array([float(x) for x in StressXX])
    StressYY = np.array([float(x) for x in StressYY])
    StressXY = np.array([float(x) for x in StressXY])
    StressTriax = np.array([float(x) for x in StressTriax])
    MisesStress = np.array([float(x) for x in MisesStress])
    ThirdINVStress = np.array([float(x) for x in ThirdINVStress])

    #Rotation angle calculation

    q = ((StressXX-StressYY)/(abs(StressXX-StressYY)))*((abs(MaxPrincStress)-abs(MinPrincStress))/(abs(abs(MaxPrincStress)-abs(MinPrincStress))))
    Beta = np.degrees(np.arctan((2*StressXY)/(StressXX-StressYY))/2)

    RotAngle = np.zeros((np.size(Beta), 1))

    for i in range(np.size(Beta)):
        if StressXX[i] == StressYY[i] and StressXY[i] != 0:
            RotAngle[i] = 45
        else:
            RotAngle[i] = 45*(1-q[i]) + q[i]*abs(Beta[i])

    RotAngle = np.array([float(x) for x in RotAngle])

    #Lode angle parameter (normalized)
    ThirdINVStress_norm = np.zeros((np.size(ThirdINVStress),1))
    for i in range(0,np.size(ThirdINVStress)):
        ThirdINVStress_norm[i] = min((ThirdINVStress[i]/MisesStress[i])**3,1)
    LodeAngle_norm = 1-(2/math.pi)*np.arccos(ThirdINVStress_norm)

    #Checking if there are points in elasticity
    ElasticPoints = np.array([x for x in range(len(EqPlastStrain)) if EqPlastStrain[x]==0])

    if ElasticPoints.size>0:
        ElastMinPrincStrain = np.array(MinPrincStrain[ElasticPoints])
        ElastMaxPrincStrain = np.array(MaxPrincStrain[ElasticPoints])
        ElastMinPrincStress = np.array(MinPrincStress[ElasticPoints])
        ElastMaxPrincStress = np.array(MaxPrincStress[ElasticPoints])
        ElasticIntPointCoordX = np.array(IntPointCoordX[ElasticPoints])
        ElasticIntPointCoordY = np.array(IntPointCoordY[ElasticPoints])
        ElastRotAngle = np.array(RotAngle[ElasticPoints])
        ElastLodeAngle_norm = np.array(LodeAngle_norm[ElasticPoints])
        ElastStressTriax = np.array(StressTriax[ElasticPoints])

        PlastMinPrincStrain = np.delete(MinPrincStrain, ElasticPoints, 0)
        PlastMaxPrincStrain = np.delete(MaxPrincStrain, ElasticPoints, 0)
        PlastMinPrincStress = np.delete(MinPrincStress, ElasticPoints, 0)
        PlastMaxPrincStress = np.delete(MaxPrincStress, ElasticPoints, 0)
        PlasticIntPointCoordX = np.delete(IntPointCoordX, ElasticPoints, 0)
        PlasticIntPointCoordY = np.delete(IntPointCoordY, ElasticPoints, 0)
        PlastRotAngle =  np.delete(RotAngle, ElasticPoints, 0)
        PlastLodeAngle_norm = np.delete(LodeAngle_norm, ElasticPoints, 0)
        PlastStressTriax = np.delete(StressTriax, ElasticPoints, 0)
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
    S_size= 10.7
    Transparency = 1
    Mark = 's'
    cjet = plt.cm.get_cmap('jet', 12)
# ------------- Eq. Plastic Strain - FEA ------------- 

    fig1 = plt.figure(figsize=(10,8))
    ax = plt.axes([0.4, 0.4, .6, .6])

    colors = EqPlastStrain
    if ElasticPoints.size>0:
        colors = np.delete(colors, ElasticPoints, 0)
        elastic = plt.scatter(ElasticIntPointCoordX,ElasticIntPointCoordY, marker=Mark, c ='grey', alpha =1, s=S_size, label=r'$\bar{\epsilon}^\mathrm{p}$')
        plastic = plt.scatter(PlasticIntPointCoordX,PlasticIntPointCoordY, marker=Mark, c =colors, cmap='jet', alpha =Transparency, s=S_size, label=r'$\bar{\epsilon}^\mathrm{p}$')
    else:
        plt.scatter(IntPointCoordX,IntPointCoordY, marker=Mark, c =colors, cmap='jet', alpha =Transparency, s=S_size, label=r'$\bar{\epsilon}^\mathrm{p}$')
                
    #plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.xlim(-30,50)
    plt.ylim(-20,80)

    
    sm1 = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=max(EqPlastStrain)))
    sm1._A = []
    tks = np.linspace(0,max(abs(EqPlastStrain)),5)
    clb1 = plt.colorbar(sm1,ticks=tks, ax=ax, shrink=0.60)
    clb1.solids.set_rasterized(True)
    clb1.set_ticklabels(['%.2f'%round(i,2) for i in tks])
    clb1.set_label(r'$\bar{\epsilon}^\mathrm{p}$',labelpad=-40,y=1.15, rotation=0)
    clb1.solids.set(alpha=1)

    plt.clim(0,max(abs(EqPlastStrain)))
    plt.axis('off')

    fig1.savefig(cwd + '\\Results\\FEA_plots\\' + test_config + '\\EqPlasticStrain_ ' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig1)


# ------------- Rotation Angle -------------

    fig2 = plt.figure(figsize=(11,8))
    ax = plt.axes([0.22, 0.24, .86, .83])
    bjet = plt.cm.get_cmap('jet', 60)
    ang = np.arange(0.5,90,1)
    
    if ElasticPoints.size>0:
        #Elastic points
        bin_elastic = np.around(ElastRotAngle, 0)

        bin_elastic_unique = Counter(bin_elastic).keys() # equals to list(set(words))
        bin_elastic_unique = np.array([float(x) for x in bin_elastic_unique])
        bin_elastic_frequency_aux= Counter(bin_elastic).values() # counts the elements' frequency
        bin_elastic_frequency_aux = np.array([float(x) for x in bin_elastic_frequency_aux])

        bin_elastic_frequency_total = np.zeros((np.size(ang), 1)) 
        for i in range(0, np.size(bin_elastic_unique)): 
            index = int(bin_elastic_unique[i])-1
            bin_elastic_frequency_total[index] = (bin_elastic_frequency_aux[i]/np.size(RotAngle))*100
        bin_elastic_frequency_total = np.array([float(x) for x in bin_elastic_frequency_total])

        #Plastic points
        bin_plastic = np.around(PlastRotAngle, 0)

        bin_plastic_unique = Counter(bin_plastic).keys() # equals to list(set(words))
        bin_plastic_unique = np.array([float(x) for x in bin_plastic_unique])
        bin_plastic_frequency_aux= Counter(bin_plastic).values() # counts the elements' frequency
        bin_plastic_frequency_aux = np.array([float(x) for x in bin_plastic_frequency_aux])
        bin_plastic_frequency = (bin_plastic_frequency_aux/np.size(RotAngle))*100

        bins_plastic_inc  = np.zeros((np.size(ang), 1))   
        bin_plastic_frequency_total = np.zeros((np.size(ang), 1)) 
        for i in range(0, np.size(bin_plastic_unique)): 
            index = int(bin_plastic_unique[i])-1
            bins_plastic_inc[index] = bin_plastic_frequency [i]/bin_plastic_frequency_aux[i]
            bin_plastic_frequency_total[index] = bin_plastic_frequency_aux[i]
        bins_plastic_inc = np.array([float(x) for x in bins_plastic_inc])

        EqPlastStrain_plastic = np.delete(EqPlastStrain, ElasticPoints, 0)
        EqPlastStrain_plastic_norm = EqPlastStrain_plastic/max(EqPlastStrain_plastic)
        colors_aux = []
        for i in range(0,np.size(bin_plastic_unique)):
            colors_aux1 = np.where(bin_plastic==bin_plastic_unique[i])
            colors_aux1 = np.array([float(x) for x in colors_aux1[0]])
            EqPlasticStrain_aux = np.zeros((np.size(colors_aux1), 1))
            for j in range(0,np.size(colors_aux1)):
                cindex = int(colors_aux1[j])
                EqPlasticStrain_aux [j] = EqPlastStrain_plastic_norm[cindex]
            colors_aux.append(max(EqPlasticStrain_aux))

        colors_inc= np.zeros((np.size(ang), 1)) 
        for i in range(0, np.size(colors_aux)): 
            index = int(bin_plastic_unique[i])-1
            colors_inc [index] = colors_aux [i]/bin_plastic_frequency_aux[i]
        colors_inc = np.array([float(x) for x in colors_inc])

        bins = np.zeros((np.size(ang), int(max(bin_plastic_frequency_aux))+1)) 
        bins_total = np.zeros((np.size(ang), int(max(bin_plastic_frequency_aux))+1))
        bins_total[:,0] = bin_elastic_frequency_total[:] 
        color_plot = np.zeros((np.size(ang), int(max(bin_plastic_frequency_aux))+1)) 
        for i in range(0,np.size(ang)):
            for k in range(0,int(max(bin_plastic_frequency_aux))):  
                if k < int(bin_plastic_frequency_total[i]):
                    bins[i,k+1] = bins_plastic_inc[i]
                    color_plot[i,k+1] = colors_inc[i]*(k+1)
                    bins_total[i,k+1] = bins_plastic_inc[i]*(k+1)+bins_total[i,0]
                if k>=int(bin_plastic_frequency_total[i]) and k!=0:
                    color_plot[i,k+1] = color_plot[i,k]
                    bins_total[i,k+1] = bins_total[i,k]

        plt.bar(ang,bin_elastic_frequency_total, width=-1, align='edge', color='grey' ,edgecolor='none',rasterized=True)
        for i in range(0,int(max(bin_plastic_frequency_aux))):
            color_set = color_plot[:,i+1]
            plt.bar(ang,bins[:,i+1], bottom=bins_total[:,i], width=-1, align='edge', color=bjet(color_set) ,edgecolor='none',rasterized=True)

    else: 
        #All points in plasticity
        bin_plastic = np.around(RotAngle, 0)

        bin_plastic_unique = Counter(bin_plastic).keys() # equals to list(set(words))
        bin_plastic_unique = np.array([float(x) for x in bin_plastic_unique])
        bin_plastic_frequency_aux= Counter(bin_plastic).values() # counts the elements' frequency
        bin_plastic_frequency_aux = np.array([float(x) for x in bin_plastic_frequency_aux])
        bin_plastic_frequency = (bin_plastic_frequency_aux/np.size(RotAngle))*100

        bins_plastic_inc  = np.zeros((np.size(ang), 1))   
        bin_plastic_frequency_total = np.zeros((np.size(ang), 1)) 
        for i in range(0, np.size(bin_plastic_unique)): 
            index = int(bin_plastic_unique[i])-1
            bins_plastic_inc[index] = bin_plastic_frequency [i]/bin_plastic_frequency_aux[i]
            bin_plastic_frequency_total[index] = bin_plastic_frequency_aux[i]
        bins_plastic_inc = np.array([float(x) for x in bins_plastic_inc])

        EqPlastStrain_norm = EqPlastStrain/max(EqPlastStrain)
        colors_aux = []
        for i in range(0,np.size(bin_plastic_unique)):
            colors_aux1 = np.where(bin_plastic==bin_plastic_unique[i])
            colors_aux1 = np.array([float(x) for x in colors_aux1[0]])
            EqPlasticStrain_aux = np.zeros((np.size(colors_aux1), 1))
            for j in range(0,np.size(colors_aux1)):
                cindex = int(colors_aux1[j])
                EqPlasticStrain_aux [j] = EqPlastStrain_norm[cindex]
            colors_aux.append(max(EqPlasticStrain_aux))

        colors_inc= np.zeros((np.size(ang), 1)) 
        for i in range(0, np.size(colors_aux)): 
            index = int(bin_plastic_unique[i])-1
            colors_inc [index] = colors_aux [i]/bin_plastic_frequency_aux[i]
        colors_inc = np.array([float(x) for x in colors_inc])

        bins = np.zeros((np.size(ang), int(max(bin_plastic_frequency_aux))+1)) 
        bins_total = np.zeros((np.size(ang), int(max(bin_plastic_frequency_aux))+1)) 
        color_plot = np.zeros((np.size(ang), int(max(bin_plastic_frequency_aux))+1)) 
        for i in range(0,np.size(ang)):
            for k in range(0,int(max(bin_plastic_frequency_aux))):
                if k < int(bin_plastic_frequency_total[i]):
                    bins[i,k+1] = bins_plastic_inc[i]
                    color_plot[i,k+1] = colors_inc[i]*(k+1)
                    bins_total[i,k+1] = bins_plastic_inc[i]*(k+1)
                if k>=int(bin_plastic_frequency_total[i]) and k!=0:
                    color_plot[i,k+1] = color_plot[i,k]
                    bins_total[i,k+1] = bins_total[i,k]

        for i in range(0,int(max(bin_plastic_frequency_aux))):
            color_set = color_plot[:,i+1]
            bars = plt.bar(ang,bins[:,i+1], bottom=bins_total[:,i], width=-1, align='edge', color=bjet(color_set) ,edgecolor='none',rasterized=True)
    
            
    sm = plt.cm.ScalarMappable(cmap=cjet,norm=plt.Normalize(vmin=0, vmax=max(EqPlastStrain)))
    sm._A = []
    tks = np.linspace(0.0,max(EqPlastStrain),9)
    clb = plt.colorbar(sm,ticks=tks, ax=ax)
    clb.solids.set_rasterized(True)
    clb.set_label(r'$\bar \epsilon ^\mathrm{P}$',labelpad=-65,y=1.07, rotation=0)
    clb.set_ticklabels(['%.2f'%round(i,2) for i in tks])
    plt.title(" ")
    plt.xlabel(r'$\gamma$ $[\degree]$' , fontsize=FS, labelpad=10)
    plt.ylabel("Density " + r'$[\%]$', fontsize=FS, labelpad=10)
    plt.xlim(0,90)
    plt.xticks(np.linspace(0,90,7))
    #plt.gca().yaxis.set_major_formatter(PercentFormatter(1, decimals=0, symbol=None))

    fig2.savefig(cwd + '\\Results\\FEA_plots\\' + test_config + '\\RotationAngle_' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig2)

# ------------- StressStates - FEA ------------- 

    fig3 = plt.figure(figsize=(10,8))
    ax = plt.axes([0.4, 0.4, .6, .6])

    colors = MinPrincStrain/MaxPrincStrain
    if ElasticPoints.size>0:
        colors = np.delete(colors, ElasticPoints, 0)
        elastic = plt.scatter(ElasticIntPointCoordX,ElasticIntPointCoordY, marker=Mark, c ='grey', alpha =1, s=S_size, label=r'$\bar{\epsilon}^\mathrm{p}$')
        plastic = plt.scatter(PlasticIntPointCoordX,PlasticIntPointCoordY, marker=Mark, c =colors, cmap='jet', alpha =Transparency, s=S_size, label=r'$\bar{\epsilon}^\mathrm{p}$')
    else:
        plt.scatter(IntPointCoordX,IntPointCoordY, marker=Mark, c =colors, cmap='jet', alpha =Transparency, s=S_size, label=r'$\bar{\epsilon}^\mathrm{p}$')
                
    #plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.xlim(-30,50)
    plt.ylim(-20,80)

    cjet = plt.cm.get_cmap('jet', 12)
    sm1 = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=min(colors),vmax=max(colors)))
    sm1._A = []
    #tks = [min(colors), -2, -1 , -0.5, max(0.6*max(colors), max(colors))] 
    tks = [min(colors), -2, -1 , -0.5, max(colors)]
    clb1 = plt.colorbar(sm1,ticks=tks, ax=ax, shrink=0.80)
    clb1.solids.set_rasterized(True)
    clb1.set_ticklabels(['%.2f'%round(i,2) for i in tks])
    clb1.set_label(r'$\epsilon_\mathrm{2}/\epsilon_\mathrm{1}$',labelpad=-40,y=1.15, rotation=0)
    clb1.solids.set(alpha=1)

    plt.clim(min(colors),max(colors))
    plt.axis('off')

    fig3.savefig(cwd + '\\Results\\FEA_plots\\' + test_config + '\\StrainStates_ ' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig3)

# ------------- Lode angle vs Stress triaxility - FEA ------------- 

    fig4 = plt.figure(figsize=(16,9))
    ax = plt.axes([0.22, 0.24, .86, .83])

    colors = EqPlastStrain
    if ElasticPoints.size>0:
        colors = np.delete(colors, ElasticPoints, 0)
        elastic = plt.scatter(ElastLodeAngle_norm,ElastStressTriax, c = 'grey', s=15,label=r'Elastic')
        plastic = plt.scatter(PlastLodeAngle_norm,PlastStressTriax, c = colors, cmap='jet', s=15,label=r'Plastic')
    else:
        plt.scatter(LodeAngle_norm,StressTriax, c = colors, cmap='jet', s=15,label=r'Plastic')

    plt.plot([-1,1],[-1/3,1/3],'k-',lw=0.8)
    x = [-1.0, 0.0, 1.0]
    y = [2.0/3.0,math.sqrt(3.0)/3.0,1.0/3.0]
    f = CubicSpline(x, y)
    xnew = np.linspace(-1, 1, num=100)
    plt.plot(xnew,f(xnew),'k-',lw=0.8)
    y = [-1.0/3.0,-math.sqrt(3.0)/3.0,-2.0/3.0]
    f = CubicSpline(x, y)
    xnew = np.linspace(-1, 1, num=100)
    plt.plot(xnew,f(xnew),'k-',lw=0.8)

    box = dict(boxstyle='Square,pad=0.05',fc='w',ec='w')
    al = 'center'
    fnt = 10

    leg = ax.legend(loc=3,frameon=False, ncol=1, fontsize=fnt, markerscale=3.0,handletextpad=0.1,bbox_to_anchor=(0.1,0.01),
                handleheight=1.5, labelspacing=0.5)

    name = 'Equibiaxial \ntension'
    plt.plot([-1,-1],[-0.75,0.75],'k-',lw=0.8,zorder=0)
    plt.text(-1,0.6,name,ha=al,va=al,ma=al,bbox=box,fontsize=fnt)

    name = 'Uniaxial \ntension'
    plt.plot([ 1, 1],[-0.75,0.75],'k-',lw=0.8,zorder=0)
    plt.text(1,0.25,name,ha=al,va=al,ma=al,bbox=box,fontsize=fnt)

    name = 'Shear'
    plt.plot([ 0, 0],[-0.75,0.75],'k--',lw=0.8,zorder=0)
    plt.text(0,0.05,name,ha=al,va=al,ma=al, bbox=box,fontsize=fnt)

    name = 'Plane strain \ncompression'
    plt.text(0,-0.51,name,ha=al,va=al,ma=al, bbox=box,fontsize=fnt)

    name = 'Plane strain \ntension'
    plt.text(0,0.51,name,ha=al,va=al,ma=al, bbox=box,fontsize=fnt)

    name = 'Uniaxial \ncompression'
    plt.text(-1,-0.25,name,ha=al,va=al,ma=al,bbox=box,fontsize=fnt)

    name = 'Equibiaxial \ncompression'
    plt.text(1,-0.6,name,ha=al,va=al,ma=al,bbox=box,fontsize=fnt)


    #plt.title(test_config + " | " + Param[i] + " sensitivity | Step:" + str(Step))
    plt.xlabel('$\\bar{\\theta}$')
    plt.ylabel('$\eta$')
    plt.ylim(-0.75,0.75)
    plt.yticks(np.linspace(-0.75,0.75,7))
    plt.xlim(-1.2,1.2)
    ax.tick_params(axis='x', colors=(0,0,0))

    cjet = plt.cm.get_cmap('jet', 12)
    sm1 = plt.cm.ScalarMappable(cmap=cjet, norm=plt.Normalize(vmin=0,vmax=max(EqPlastStrain)))
    sm1._A = []
    tks = np.linspace(0,max(abs(EqPlastStrain)),5)
    clb1 = plt.colorbar(sm1,ticks=tks, ax=ax)
    clb1.solids.set_rasterized(True)
    clb1.set_ticklabels(['%.2f'%round(i,2) for i in tks])
    clb1.set_label(r'$\bar{\epsilon}^\mathrm{p}$',labelpad=-40,y=1.1, rotation=0)
    clb1.solids.set(alpha=1)

    plt.clim(0,max(abs(EqPlastStrain)))

    fig4.savefig(cwd + '\\Results\\FEA_plots\\' + test_config + '\\LodeAngleStressTriax_ ' + test_config + '.jpg', dpi=600, bbox_inches='tight', pad_inches=.15)
    plt.close(fig4)

    if show_plots == 'yes':
        plt.show()
    else:
        pass
        
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()