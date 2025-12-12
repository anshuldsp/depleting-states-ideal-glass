# Module importieren

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import seaborn as sns
import pandas as pd
import math

from csv import writer
from typing import List
from scipy.stats import linregress
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
mpl.use('TkAgg')

import scienceplots

#mpl.rcParams["axes.spines.right"] = False
#mpl.rcParams["axes.spines.top"] = False

# Eigene Module

import Distribution as dst

# Class Definition

class E_IS:
    '''Reads an energy distribution G(E) and calculates different properties.'''
    def __init__(self,dist:dst.Distribution):
        self.distribution = dist
        self.sigma = dist.var
        dist.calc_all()
        self.ener = dist.scale
        self.ener_dist = dist.value

        self.temp = []
        self.avg = []
        self.rate = []
        self.app = []
        self.target = 0.0
        return
    
    def set_temp(self,temp_min:float=0.05,temp_max:float=10.0,temp_step:float=0.0005) -> None:
        '''Sets the investigated temperature interval.\\
        For binomial distributions: temp_min = 0.006, temp_max = 10.0, temp_step = 0.0005\\
        For gaussian distributions: temp_min = 0.06, temp_max = 10.0, temp_step = 0.005'''
        self.temp = np.arange(temp_min,temp_max,temp_step)
        for temp in self.temp:
            if temp <= 0.0:
                raise ValueError("Negative or 0K")
        return
    
    def calc_avg(self,temp:float=0.1) -> float:
        '''Calculates the average energy <E> for a given temperature.'''
        inv_temp = 1/temp
        ener_avg = 0.0
        part_sum = 0.0
        for i,e_i in enumerate(self.ener):
            exponent = np.log(self.ener_dist[i])-e_i*inv_temp
            ener_avg += e_i*np.exp(exponent)
            part_sum += np.exp(exponent)
        return (ener_avg/part_sum)
    
    def calc_rate(self,temp:float=0.1) -> float:
        '''Calculates the rate \Gamma for a given temperature.'''
        inv_temp = 1/temp
        rate = 0.0
        part_sum = 0.0
        for i,e_i in enumerate(self.ener):
            exponent = np.log(self.ener_dist[i])-e_i*inv_temp
            rate += self.ener_dist[i]
            part_sum += np.exp(exponent)
        return np.log(rate/part_sum)
    
    def calc_app(self) -> None:
        '''Calculates the apparent activation energy for a given temperature.\\
        NOTE: It can be shown that this is equal to -<E>.'''
        if len(self.rate) != len(self.temp):
            self.rate = [self.calc_rate(temp) for temp in self.temp]
        self.app = -1.0*np.gradient(self.rate,[(1/temp) for temp in self.temp])
        return
    
    def calc_target(self) -> None:
        '''Calculates the onset temperature.'''
        if len(self.avg) != len(self.temp):
            self.avg = [self.calc_avg(temp) for temp in self.temp]
        x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.avg[i])]
        y = [avg for avg in self.avg if np.isfinite(avg)]
        deriv = np.gradient(np.gradient(y,x),x)

        max = 0.0
        for i,element in enumerate(deriv):
            if element >= max:
                max = element
        check = 2.0*max/3.0
        counter = len(deriv)-1
        for element in reversed(deriv):
            if element >= check:
                target = counter
                break
            counter -= 1
        self.target = x[target]
        return

    def calc_all(self) -> None:
        '''Gathers all necessary calc_*-functions.'''
        self.set_temp()
        self.calc_target()
        self.calc_app()
        return
    
    def plot_paper(self) -> None:
        '''Plots the figure 5 of the PRX Paper.'''
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'sans-serif'
        plt.style.use(['science','nature'])

        target = 5

        fig = plt.figure(figsize=(3.39,2.0))

        self.calc_all()

        x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.avg[i])]
        x = x[::2]

        y1 = [app for app in self.app if np.isfinite(app)]
        y1 = y1[::2]

        sub1 = fig.add_subplot(122)
        sub1.set_ylim(-0.5,12)
        sub1.set_xlabel(r"$1/T$")
        sub1.plot(x[::5],y1[::5],marker='o',markersize=2,linestyle='None')
        
        sel_temp = x[-10:]
        sel_val = y1[-10:]
        slope,intercept,r_value,p_value,std_err = linregress(sel_temp,sel_val)
        sub1.plot(x,[slope*i+intercept for i in x],linestyle="--",color="black")

        sub1.spines['right'].set_visible(False)
        sub1.spines['top'].set_visible(False)
        sub1.xaxis.set_ticks_position('bottom')
        sub1.yaxis.set_ticks_position('left')
        sub1.axvline(target,linestyle="--",color="red")

        ins1 = inset_axes(sub1,width="40%",height="40%",loc="lower right",bbox_to_anchor=(0,0.185,1,1),bbox_transform=sub1.transAxes,borderpad=0)
        ins1.set_title(r"$-\ln\langle\Gamma\rangle$",fontsize=7.0)
        ins1.set_xlabel(r"$1/T$",labelpad=0)
        y2 = [(-i) for i in self.rate]
        y2 = y2[::2]
        ins1.plot(x,y2,marker="o",markersize=1,linestyle="None")
        ins1.set_xlim(0,20)
        ins1.axvline(target,linestyle="--",color="red")

        slope2,intercept2,r_value,p_value,std_err = linregress(x[44:47],y2[44:47])
        ins1.plot(x,[slope2*x+intercept2 for x in x],linestyle="--",color="black")

        #----------------

        y3 = [-1.0*i for i in y1]

        sub2 = fig.add_subplot(121)
        sub2.plot(x,y3,marker="o",markersize=2,linestyle="None")
        sub2.plot(x,[-(slope*i+intercept) for i in x],linestyle="--",color="black")
        sub2.set_ylim(-12,0.5)
        sub2.set_xlabel(r"$1/T$")

        sub2.spines['right'].set_visible(False)
        sub2.spines['top'].set_visible(False)
        sub2.xaxis.set_ticks_position('bottom')
        sub2.yaxis.set_ticks_position('left')
        sub2.axvline(target,linestyle="--",color="red")

        y4 = [(avg + (1/self.temp[i])*self.sigma*self.sigma) for i,avg in enumerate(self.avg) if np.isfinite(self.avg[i])]
        y4 = y4[::2]

        ins2 = inset_axes(sub2,width="40%",height="40%",loc="center right",bbox_to_anchor=(0,0.15,1,1),bbox_transform=sub2.transAxes,borderpad=0)
        ins2.set_title(r"$\langle E_{IS}\rangle-E_{IS,Gauss}$",fontsize=7.0)
        ins2.set_ylim(-0.5,8)
        ins2.set_xlim(0,20)
        ins2.plot(x,y4,marker="o",markersize=1,linestyle="None")
        ins2.set_xlabel(r"$1/T$",labelpad=0)
        ins2.axhline(0,linestyle="--",color="black")
        ins2.axvline(target,linestyle="--",color="red")

        sub1.text(-5.66, 11.785, r'${\bf b.}$')
        sub2.text(-5, 0.215, r'${\bf a.}$')

        sub2.text(1, 0.215, r'$\langle E_{IS} \rangle$')
        sub1.text(1.75, 11.785, r'$E_{app}$')

        print(plt.rcParams["font.size"])

        plt.tight_layout()
        plt.savefig("dist_150.pdf",dpi=400)
        plt.savefig("dist_150.png")
        plt.show()
        plt.close()      
        return

        


    def plot_prop(self,mode:int=0,target:bool=False,taylor:bool=False,path:str="test") -> None:
        '''Plots the calculated properties.\\
        mode = 0: Average Energy.\\
        mode = 1: Rate.\\
        mode = 2: Apparent activation energy.
        mode = 3: <E> - <E> (Gauss).'''

        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'sans-serif'

        plt.style.use(['science','nature'])


        fig = plt.figure(figsize=(3.39,2.0))


        self.calc_all()
        if mode == 0:
            x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.avg[i])]
            y = [avg for avg in self.avg if np.isfinite(avg)]
            path = "avg"
        elif mode == 1:
            x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.rate[i])]
            y = [rate for rate in self.rate if np.isfinite(rate)]
            path = "rate"
        elif mode == 2:
            x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.app[i])]
            y1 = [app for app in self.app if np.isfinite(app)]
            y2 = [(avg + (1/self.temp[i])*self.sigma*self.sigma) for i,avg in enumerate(self.avg) if np.isfinite(self.avg[i])]
        
            path = "app"
        elif mode == 3:
            x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.avg[i])]
            y = [(avg + (1/self.temp[i])*self.sigma*self.sigma) for i,avg in enumerate(self.avg) if np.isfinite(self.avg[i])]
        
        ax = fig.add_subplot(122)
        ax.set_ylim(-0.5,12)

        p1 = ax.plot(x[::5],y1[::5],marker='o',markersize=2,linestyle="None")
        selected_temp = x[-10:]
        selected_value = y1[-10:]
        slope,intercept,r_value,p_value,std_err = linregress(selected_temp,selected_value)
            
        x = x[::2]
        y1 = y1[::2]

        p2 = ax.plot(x,[slope*i+intercept for i in x],linestyle="--",color="black")#color="#1f77b4")

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        

        
        ax.set_xlabel(r"$1/T$")

        axins2 = inset_axes(ax,width="40%",height="40%",
                            loc="lower right",
                            bbox_to_anchor=(0,0.15,1,1),
                            bbox_transform=ax.transAxes,
                            borderpad=0)
        axins2.set_title(r"$\ln\langle\Gamma\rangle(1/T)$")
        axins2.set_xlabel(r"(1/T)")
        y4 = [(i) for i in self.rate]
        y4 = y4[::2]
        axins2.plot(x,y4,marker="o",markersize=2,linestyle="None")
        axins2.set_xlim(0,20)

        slope2,intercept2,r_value,p_value,std_err = linregress(x[88:94],y4[88:94])
        axins2.plot(x, [slope2*x + intercept2 for x in x], 
                      color="black", linestyle="--")
        
        print(slope2)



        y3 = [-1.0*i for i in y1]
        
        if mode == 0:
            ax.set_ylabel(r"$\langle E\rangle$")
            ax.set_xlim(0,50)
            ymin = y1[0]
            ymax = y1[-1]
            ax.set_ylim(ymin,ymax)
        elif mode == 1:
            ax.set_ylabel(r"$\langle \Gamma \rangle$")
            ax.set_xlim(0,self.distribution.M)
            ymin = -self.distribution.M*0.5
            ymax = 0.0
            ax.set_ylim(ymin,ymax)
        elif mode == 2: 
            ax.set_ylabel(r"$E_{app}$")
            ax.set_xlim(0,20)
            ymin = y1[-1]
            ymax = y1[0]
            ax.set_ylim(ymin,ymax)
        elif mode == 3:
            ax.set_ylabel(r"$e_{IS}-E_{IS,Gauss}$",fontsize=15)
            ax.axhline(0,color="red",linestyle="--")
            ymin = y1[-1]
            ymax = y1[0]

        
        ax = fig.add_subplot(121)

        axins = inset_axes(ax,width="40%",height="40%",loc="upper right")
        ax.plot(x,y3,marker='o',markersize=2,linestyle='None')
        ax.plot(x,[-(slope*i+intercept) for i in x],linestyle="--",color="black")
        axins.axvline(2.5,linestyle="--",color="red")
        ax.set_ylim(-12,0.5)

        ax.text(14, 0.215, r'$\bf{b.}$')

        axins.set_title(r"$\langle E_{IS}\rangle-E_{IS,Gauss}$")

        p1 = axins.plot(x,y2,marker='o',markersize=1,linestyle='None')
        axins.axhline(0,linestyle="--",color="black")
        ax.set_xlabel(r"$1/T$")
        axins.set_ylim(-0.5,12)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.axvline(2.5,linestyle="--",color="red")
        plt.savefig("dist_150.pdf")
        plt.savefig("dist_150.png")
        plt.show()
        plt.close()      
        return
    
    def write_prop(self) -> None:
        '''Writes the calculated properties.\\
        NOTE: The name of the output file is hardcoded.'''
        self.calc_all()
        path = f"eis_{self.distribution.M}.csv"
        with open(path,"w",newline="") as csvfile:
            csv_writer = writer(csvfile)
            csv_writer.writerow([f"1/T = {self.target}","avg","rate","app"])
            for i in range(len(self.temp)):
                csv_writer.writerow([
                    str(1/self.temp[i]),
                    str(self.avg[i]),
                    str(self.rate[i]),
                    str(self.app[i])
                ])
        csvfile.close()
        return

    
def Test(M:int=50,var:float=1.0):
    dist0 = dst.Distribution(M=M,var=var)
    #dist0.plot_dist()
    eis0 = E_IS(dist0)
    eis0.plot_paper()
    #eis0.plot_prop(mode=3,target=True)
    #eis0.plot_prop(mode=0,target=True,taylor=True)
    #eis0.plot_prop(mode=0,target=False,taylor=True)
    #eis0.plot_prop(mode=1,target=True)
    #eis0.plot_prop(mode=2,target=True)
    #eis0.write_prop()

#Test(M=3)
#Test(M=10)
Test(M=150)
#Test(M=150)
#Test(M=200)
