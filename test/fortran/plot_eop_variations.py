#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np

vmjd=[]; vxpot=[]; vypot=[]; vdut1ot=[]; vxplb=[]; vyplb=[]; vdut1lb=[];

if __name__ == "__main__":
    for line in sys.stdin:
        if line[0] != '#':
            print(line)
            mjd, xpot, ypot, dut1ot, xplb, yplb, dut1lb, dlodlb = [ float(x) for x in line.split() ]
            vmjd.append(mjd); 
            vxpot.append(xpot);vypot.append(ypot);vdut1ot.append(dut1ot); 
            vxplb.append(xplb);vyplb.append(yplb);vdut1lb.append(dut1lb); 

    ## do nothing on empty input
    if len(vmjd) <= 1: sys.exit(1)

    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(vmjd, vxpot, vmjd, vxplb)
    axs[0, 0].set_title(r'$\Delta x_p$')
    axs[0, 0].legend([r'$\Delta x_p$ ocean tide', r'$\Delta x_p$ libration'])
    
    axs[0, 1].plot(vmjd, vypot, vmjd, vyplb)
    axs[0, 1].set_title(r'$\Delta y_p$')
    axs[0, 1].legend([r'$\Delta y_p$ ocean tide', r'$\Delta y_p$ libration'])
    
    axs[1, 0].plot(vmjd, vdut1ot, vmjd, vdut1lb)
    axs[1, 0].set_title(r'$\Delta UT1$')
    axs[1, 0].legend([r'$\Delta UT1$ ocean tide', r'$\Delta UT1$ libration'])

    plt.show()


