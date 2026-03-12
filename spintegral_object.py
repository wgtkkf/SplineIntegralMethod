# Third-order spline interporation & integral calculation
# Change input files name & input files number
# Read gap & flux data files
# Coded by Takuro TOKUNAGA
# Last modified: May 15 2018
# Updated: March 10, 11, 2026

import numpy as np
import time
import pandas as pd
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

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
        sa = self.f1           # small a
        sh = self.s2 - self.s1 # small h
        sb = (1/sh)*(self.f2-self.f1) - (sh/3)*(2*self.v1 + self.v2) # small b
        sc = self.v1                    # small c
        sd = (self.v2 - self.v1)/(3*sh) # small d

        # define of interporation function
        ipfx = sa + sb*(x-self.s1) + sc*np.power((x-self.s1),2) + sd*np.power((x-self.s1),3)

        return ipfx        

class SplineIntegral:
    def __init__(self, ):


    def loop_calculation(self, ):


# parameters
number = 3700 # number of data points x & fx in flux*.txt, change here, for Gold
filenum = 10 # number of flux txt files, not the inside of txt file
integral = 0 # initialization
sum = 0      # initialization

#
sx=np.zeros((number), dtype='float64')
fx=np.zeros((number), dtype='float64')

#  matrix and vectors
ma=np.zeros((number,number), dtype='float64')
vx=np.zeros((number), dtype='float64')
vb=np.zeros((number), dtype='float64')

# vector for gap
vg=np.zeros((filenum), dtype='float64') # vector for gap


# output file
f = open('../fe/power_test.txt', 'w')

def main():
    # main start
    start = TimeCounter(time.time())
    msg = Comments('Calculation started.', 'Calculation completed.')

    # message method
    msg.begin()

    # file read x & fx
    for j in range(0, filenum):    
        data = pd.read_csv("../fe/flux/au/flux"+str(j)+".txt", sep=" ", header=None) # change file name here tab:\t for CHPC
        data.columns = ["omega", "prop", "evan", "prop+evan"]

        gdata = pd.read_csv("../fe/flux/spline/gap.txt", sep=" ", header=None) # gap data
        gdata.columns = ["gap"]
        vg[j] = gdata.iat[j,0] # 0th line

        # input data into tables
        for i in range(0, number):
            sx[i] = data.iat[i,0] # x line
            fx[i] = data.iat[i,3] # fx line, 1:prop, 2:evan, 3:total

        # initialization of matrix & vectors
        # matrix A
        for i in range(0, number): # from 0 to 13, total 14
            # sub-diagonal components
            if i>0 and i<number-1: # 1 to 12, number=14
                ma[i][i-1] = sx[i] - sx[i-1] # left
                ma[i][i+1] = sx[i+1] - sx[i] # right

            #diagonal components
            if i==0:
                ma[i][i] = 1
            elif i>0 and i<number-1: # 1 to 12
                ma[i][i] = 2*(ma[i][i-1]+ma[i][i+1])
            elif i==number-1: # i=13
                ma[i][i] = 1

        # vector b
        for i in range(0, number): # from 0 to 13, total 14
            if i>0 and i<number-1: # from 1 to 12
                h1 = sx[i+1]-sx[i]
                h0 = sx[i]-sx[i-1]
                a2 = fx[i+1]
                a1 = fx[i]
                a0 = fx[i-1]

                vb[i] = (3/h1)*(a2-a1)-(3/h0)*(a1-a0)

            # reset parameters
            h1 = 0
            h0 = 0
            a2 = 0
            a1 = 0
            a0 = 0        

        # calculation of vector x, LU method
        vx = np.linalg.solve(ma, vb) # maybe this is correct

        # reset parameters for next files
        sum = 0

        # integral calculation
        for i in range(0, number-1):

            # x, fx(i), fx(i+1), sx(i), sx(i+1), vx(i), vx(i+1)
            ipf = InterPolation(fx[i], fx[i+1], sx[i], sx[i+1], vx[i], vx[i+1])

            integral = quad(ipf.function, sx[i], sx[i+1])[0]
            sum += integral

        # output for file
        f.write(str(vg[j]))
        f.write(str(' '))
        f.write(str(sum)) # [W/m2]
        f.write('\n')

    # file close
    f.close()

    # message method
    msg.end()

    # time display
    end = TimeCounter(time.time())
    print(f"### elapsed_time: {end.time_return() - start.time_return():.2f} [sec] ###")

# --- main routine ---
if __name__ == '__main__':
    main()
