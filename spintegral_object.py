## About this code
# Third-order spline interporation & integral calculation
# Read files of heat flux file and gap externally
# Coded by Takuro Tokunaga

## Update
# I updated my scripting style in 2018 to an object-oriented style in 2026

## Update history
# Last modified: May 15 2018
# Updated: March 10-12, 2026

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.integrate import trapz, simps, quad, quadrature, romberg
import time

class Comments:
    def __init__(self, comment1, comment2):
        self.com1 = comment1
        self.com2 = comment2

    def begin(self):
        print ("### " + str(self.com1)  + "     ###")

    def end(self):
        print ("### " + str(self.com2)  + "   ###")

    pass

class TimeCounter:
    def __init__(self, time_counter):
        self.time_counter = time_counter        
        return None

    def time_return(self):
        return self.time_counter
    
## This class was generated supported by Gemini Pro ##
class TextFileAnalyzer:
    def __init__(self, folder_path: str, folder_path_gap: str):
        # Convert the string path to a Path object when the class is initialized
        self.folder = Path(folder_path)
        self.file_path = Path(folder_path_gap)

    def count_txt_files(self) -> int:
        """Counts how many .txt files are in the specified folder."""
        if not self.folder.is_dir():
            print(f"Error: The directory '{self.folder}' does not exist.")
            return 0
        
        # .glob("*.txt") finds all files ending in .txt
        txt_files = list(self.folder.glob("*.txt"))
        return len(txt_files)

    def count_lines_in_file(self) -> int:
        """Counts the total number of lines in the initialized file."""
        # No need to combine folder + filename anymore!
        if not self.file_path.is_file():
            print(f"Error: The file '{self.file_path}' does not exist.")
            return 0
            
        try:
            with self.file_path.open('r', encoding='utf-8') as file:
                return sum(1 for _ in file)
        except UnicodeDecodeError:
            print(f"Error: Could not read '{self.file_path.name}'. It might not be a standard text file.")
            return 0

class InterPolation:
    def __init__(self, f1, f2, s1, s2, v1, v2): # fx(i), fx(i+1), sx(i), sx(i+1), vx(i), vx(i+1)        
        self.f1 = f1
        self.f2 = f2
        self.s1 = s1
        self.s2 = s2
        self.v1 = v1
        self.v2 = v2
        return None

    def function(self, x):
        # coefficients
        sa = self.f1                                                 # small a
        sh = self.s2 - self.s1                                       # small h
        sb = (1/sh)*(self.f2-self.f1) - (sh/3)*(2*self.v1 + self.v2) # small b
        sc = self.v1                                                 # small c
        sd = (self.v2 - self.v1)/(3*sh)                              # small d

        # define of interporation function
        ipfx = sa + sb*(x-self.s1) + sc*np.power((x-self.s1),2) + sd*np.power((x-self.s1),3)

        return ipfx

class SplineIntegral:
    def __init__(self, dpoint, filenum, file_path_gap: str, file_path_flux: str):

        # paths
        self.file_path_gap = file_path_gap
        self.file_path_flux = file_path_flux

        # parameters
        self.dpoint = dpoint   # dpoint of data points x & fx in flux*.txt, change here, for Gold

        self.filenum = filenum # dpoint of flux txt files, not the inside of txt file
        self.integral = 0      # initialization
        self.sum = 0           # initialization

        #
        self.sx=np.zeros((self.dpoint), dtype='float64')
        self.fx=np.zeros((self.dpoint), dtype='float64')

        # vector for gap
        self.vg=np.zeros((self.filenum), dtype='float64') # vector for gap

    def matrix_operation(self):

        #  matrix and vectors
        self.ma=np.zeros((self.dpoint, self.dpoint), dtype='float64')
        self.vx=np.zeros((self.dpoint), dtype='float64')
        self.vb=np.zeros((self.dpoint), dtype='float64')        

        # matrix A
        for i in range(0, self.dpoint):
            # sub-diagonal components
            if i > 0 and i < self.dpoint - 1:
                self.ma[i][i-1] = self.sx[i] - self.sx[i-1] # left
                self.ma[i][i+1] = self.sx[i+1] - self.sx[i] # right

            #diagonal components
            if i==0:
                self.ma[i][i] = 1
            elif i > 0 and i < self.dpoint-1:
                self.ma[i][i] = 2*(self.ma[i][i-1] + self.ma[i][i+1])
            elif i==self.dpoint - 1:
                self.ma[i][i] = 1

        # vector b
        for i in range(0, self.dpoint):
            if i > 0 and i < self.dpoint - 1:
                h1 = self.sx[i+1] - self.sx[i]
                h0 = self.sx[i] - self.sx[i-1] 
                a2 = self.fx[i+1]
                a1 = self.fx[i]
                a0 = self.fx[i-1]

                self.vb[i] = (3/h1)*(a2-a1)-(3/h0)*(a1-a0)                    

        # calculation of vector x, LU method
        vx = np.linalg.solve(self.ma, self.vb)

        return vx

    def loop_calculation(self):
        # output file
        f = open('../fe/power.txt', 'w')

        # file read x & fx
        for j in range(0, self.filenum):    
            gdata = pd.read_csv(self.file_path_gap, sep=" ", header=None) # gap data
            gdata.columns = ["gap"]
            self.vg[j] = gdata.iat[j,0] # 0th line

            data = pd.read_csv(self.file_path_flux + str("/flux") + str(j) + ".txt", sep=" ", header=None) # change file name here, add tab:\t for CHPC
            data.columns = ["omega", "prop", "evan", "prop + evan"]            

            # input data into tables
            for i in range(0, self.dpoint):                
                # matrix operation
                self.sx[i] = data.iat[i, 0]
                self.fx[i] = data.iat[i, 3]            
            
            self.matrix_operation() # x line, fx line, 1:prop, 2:evan, 3:total

            # integral calculation
            for i in range(0, self.dpoint - 1):

                # x, fx(i), fx(i+1), sx(i), sx(i+1), vx(i), vx(i+1)
                integral = quad(InterPolation(self.fx[i], self.fx[i+1], self.sx[i], self.sx[i+1], self.vx[i], self.vx[i+1]).function, self.sx[i], self.sx[i+1])[0]
                self.sum += integral

            # output for file
            f.write(str(self.vg[j]))
            f.write(str(' '))
            f.write(str(self.sum)) # [W/m2]
            f.write('\n')

            # reset parameters for next files
            self.sum = 0

        # file close
        f.close()

def main():
    # main start

    ## path input
    flux_file_path = "../fe/flux/au/flux0.txt"
    flux_path = "../fe/flux/au"
    gap_file_path = "../fe/flux/spline/gap.txt"    
    ##

    start = TimeCounter(time.time())
    msg = Comments('Calculation started.', 'Calculation completed.')    

    analyzer = TextFileAnalyzer(flux_path, flux_file_path)
    total_txt_files = analyzer.count_txt_files()    
    total_lines = analyzer.count_lines_in_file()

    # message method
    msg.begin()

    # Spline integral class
    si = SplineIntegral(total_lines, total_txt_files, gap_file_path, flux_path) # si: Spline Integral    
    si.loop_calculation()

    # message method
    msg.end()

    # time display
    end = TimeCounter(time.time())
    print(f"### elapsed_time: {end.time_return() - start.time_return():.2f} [sec] ###")

# --- main routine ---
if __name__ == '__main__':
    main()
