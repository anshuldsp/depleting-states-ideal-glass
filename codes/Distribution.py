# Import Modules

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
from csv import writer
from typing import List

# Class Definitions

class Distribution:
    '''Creates a binomial distribution with given M (M over i) and variance (var).'''
    def __init__(self,M:int=3,var:float=1.0) -> None:
        '''Initializes the binomial distribution with given M (M over i) and variance (var).'''
        self.M = M
        self.var = var
        self.bins = [i for i in range(M+1)] 
        self.scale = []
        self.value = []
        return
    
    def calc_scale(self) -> None:
        '''Calculates the scale to have a mean of 0.0 and the given variance.'''
        if len(self.scale) != len(self.bins):
            for bin in self.bins:
                c = (2*self.var)/(float(self.M)**0.5)
                scale = (bin-0.5*float(self.M))*c
                #scale = (bin-(0.5*float(self.M)-200))*c
                self.scale.append(scale)
            return
        else:
            return
        
    def calc_value(self) -> None:
        '''Calculates the value of the binomial distribution.'''
        if len(self.value) != len(self.bins):
            for bin in self.bins:
                a = math.factorial(self.M)/(math.factorial(bin)*math.factorial(self.M-bin))
                a = a*2**float(-self.M)
                self.value.append(a)
            return
        else:
            return
        
    def calc_all(self) -> None:
        '''Gathers all necessary calc_*-Functions.'''
        self.calc_value()
        self.calc_scale()

    def plot_dist(self) -> None:
        '''Plots the binomial distribution. \\
        NOTE: The name of the output file is hardcoded.'''
        self.calc_all()
        y = self.value
        x = self.scale

        sns.set_context("paper",font_scale=1.0,rc={"line.linewidth":4})
        plt.figure()

        ax = sns.lineplot(x=x,y=y,color="#140B34",
                          linestyle="--",markers=True,marker="o")
        ax.grid(True,which="major",linestyle="-",linewidth=1.3,color="lightgray")
        ax.label_outer()
        ax.set_xlabel("E")
        ax.set_ylabel("G(E)")
        ax.set_title(f"M = {self.M}")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
            spine.set_color("black")
        ax.set_axisbelow(True)
        ax.minorticks_on()

        plt.tight_layout()
        plt.savefig(f"dist_{self.M}.pdf")
        plt.savefig(f"dist_{self.M}.png")
        plt.show()
        plt.close()
        
        return
    
    def plot_comp(self,other) -> None:
        '''Plots two binomial distributions for comparison.\\
        NOTE: The name of the output file is hardcoded.'''
        self.calc_all()
        other.calc_all()
        y0 = self.value
        x0 = self.scale
        y1 = other.value
        x1 = other.scale

        sns.set_context("paper",font_scale=1.0,rc={"line.linewidth":4})
        plt.figure()

        ax = sns.lineplot(x=x0,y=y0,color="#440154FF",
                          linestyle="--",markers=True,marker="o",
                          label=f"M = {self.M}")
        ax = sns.lineplot(x=x1,y=y1,color="#95D840FF",
                          linestyle="--",markers=True,marker="o",
                          label=f"M = {other.M}")
        ax.grid(True,which="major",linestyle="-",linewidth=1.3,color="lightgray")
        ax.label_outer()
        ax.set_xlabel("E")
        ax.set_ylabel("G(E)")
        ax.set_title(f"M = {self.M}")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
            spine.set_color("black")
        ax.set_axisbelow(True)
        ax.minorticks_on()
        ax.legend()

        plt.tight_layout()
        plt.savefig(f"comp_{self.M}_{other.M}.pdf")
        plt.savefig(f"comp_{self.M}_{other.M}.png")
        plt.show()
        plt.close()

        return
    
    def write_dist(self) -> None:
        '''Writes the binomial distribution. \\
        NOTE: The name of the output file is hardcoded.'''
        self.calc_all()
        path = f"dist_{self.M}.csv"
        with open(path,"w",newline="") as csvfile:
            csv_writer = writer(csvfile)
            csv_writer.writerow(["i","E","G(E)"])
            for i in range(len(self.bins)):
                csv_writer.writerow([
                    str(self.bins[i]),
                    str(self.scale[i]),
                    str(self.value[i])
                ])
        csvfile.close()
        return
    
def Test(M:int=3,var:float=1.0):
    '''Tests the binomial distribution class.'''
    dist0 = Distribution(M=M,var=var)
    dist0.plot_dist()

    dist1 = Distribution(M=50,var=var)
    dist1.plot_comp(dist0)
