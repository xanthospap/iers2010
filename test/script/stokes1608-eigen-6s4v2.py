#! /usr/bin/python

import math
import datetime
import sys

## ---------------------------------------------------------------------------
## Dummy script to compute C(3,2) and S(3,2) for given epochs using the 
## EIGEN-6S4v2 model.
## This script is meant to assist testing of parsing ICGEM files.
## ---------------------------------------------------------------------------

def d2s(date_str):
    """ handle dates with minutes == 60, e.g. 20041226.0060 """
    date_str = date_str.strip()
    try:
        return datetime.datetime.strptime(date_str, '%Y%m%d.%H%M')
    except:
        pass
    if date_str[-2:] == '60':
        return datetime.datetime.strptime(date_str[:-2], '%Y%m%d.%H') + datetime.timedelta(hours=1)
    else:
        print('ERROR parsing date: {:}'.format(date_str), file=sys.stderr)
        sys.exit(1)
        

# Model EIGEN-6S4v2 for (n,m) = (16,8)
eigen_6s4 = [
{'coef':'gfct', 'C': -2.12058698852E-08, 'S':  5.40872023682E-09, 't0':  '19500101.0000', 't1': '20020815.0817'},
{'coef':'trnd', 'C':  0.00000000000E+00, 'S':  0.00000000000E+00, 't0':  '19500101.0000', 't1': '20020815.0817'},
{'coef':'acos', 'C':  3.04472200919E-13, 'S':  4.34034097939E-12, 't0':  '19500101.0000', 't1': '20020815.0817', 'per': 1.0},
{'coef':'asin', 'C': -1.53913671236E-12, 'S':  2.57940752999E-12, 't0':  '19500101.0000', 't1': '20020815.0817', 'per': 1.0},
{'coef':'acos', 'C': -1.06616774764E-12, 'S': -5.25926073450E-14, 't0':  '19500101.0000', 't1': '20020815.0817', 'per': 0.5},
{'coef':'asin', 'C':  1.33608393287E-12, 'S': -1.46418050564E-13, 't0':  '19500101.0000', 't1': '20020815.0817', 'per': 0.5},
{'coef':'gfct', 'C': -2.12058698852E-08, 'S':  5.40872023682E-09, 't0':  '20020815.0817', 't1': '20030101.0000'},
{'coef':'trnd', 'C':  3.79790903585E-12, 'S':  2.02242793574E-12, 't0':  '20020815.0817', 't1': '20030101.0000'},
{'coef':'acos', 'C':  1.60851463720E-12, 'S':  5.03613224519E-12, 't0':  '20020815.0817', 't1': '20030101.0000', 'per': 1.0},
{'coef':'asin', 'C': -1.76110217948E-12, 'S':  5.08228762235E-12, 't0':  '20020815.0817', 't1': '20030101.0000', 'per': 1.0},
{'coef':'acos', 'C': -3.26111296453E-12, 'S':  1.17049226939E-13, 't0':  '20020815.0817', 't1': '20030101.0000', 'per': 0.5},
{'coef':'asin', 'C':  2.26953098781E-12, 'S':  2.06352311836E-12, 't0':  '20020815.0817', 't1': '20030101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12044271476E-08, 'S':  5.40948851026E-09, 't0':  '20030101.0000', 't1': '20040101.0000'},
{'coef':'trnd', 'C': -1.06666501433E-12, 'S': -3.06700848133E-12, 't0':  '20030101.0000', 't1': '20040101.0000'},
{'coef':'acos', 'C':  4.57634266678E-13, 'S':  4.96570369982E-12, 't0':  '20030101.0000', 't1': '20040101.0000', 'per': 1.0},
{'coef':'asin', 'C': -1.01922675513E-13, 'S':  2.89465354803E-12, 't0':  '20030101.0000', 't1': '20040101.0000', 'per': 1.0},
{'coef':'acos', 'C': -2.09021164630E-12, 'S':  2.30126849951E-13, 't0':  '20030101.0000', 't1': '20040101.0000', 'per': 0.5},
{'coef':'asin', 'C':  1.55630382161E-12, 'S':  1.32749693272E-13, 't0':  '20030101.0000', 't1': '20040101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12054938127E-08, 'S':  5.40642150178E-09, 't0':  '20040101.0000', 't1': '20041226.0060'},
{'coef':'trnd', 'C':  4.68357414841E-13, 'S':  2.53014156071E-12, 't0':  '20040101.0000', 't1': '20041226.0060'},
{'coef':'acos', 'C': -6.25842073093E-13, 'S':  4.56159324176E-12, 't0':  '20040101.0000', 't1': '20041226.0060', 'per': 1.0},
{'coef':'asin', 'C': -1.83832940575E-12, 'S':  2.80215749805E-12, 't0':  '20040101.0000', 't1': '20041226.0060', 'per': 1.0},
{'coef':'acos', 'C': -9.14771633965E-13, 'S':  7.93545034377E-13, 't0':  '20040101.0000', 't1': '20041226.0060', 'per': 0.5},
{'coef':'asin', 'C':  1.50568221945E-12, 'S':  1.27627716621E-12, 't0':  '20040101.0000', 't1': '20041226.0060', 'per': 0.5},
{'coef':'gfct', 'C': -2.12038045127E-08, 'S':  5.40964566611E-09, 't0':  '20041226.0060', 't1': '20060101.0000'},
{'coef':'trnd', 'C': -2.20988010217E-12, 'S': -4.25201313231E-12, 't0':  '20041226.0060', 't1': '20060101.0000'},
{'coef':'acos', 'C':  7.59581151812E-14, 'S':  4.35521073170E-12, 't0':  '20041226.0060', 't1': '20060101.0000', 'per': 1.0},
{'coef':'asin', 'C': -1.58101921481E-12, 'S':  4.97767473143E-13, 't0':  '20041226.0060', 't1': '20060101.0000', 'per': 1.0},
{'coef':'acos', 'C': -7.97308643080E-13, 'S':  6.49151243715E-15, 't0':  '20041226.0060', 't1': '20060101.0000', 'per': 0.5},
{'coef':'asin', 'C':  1.69014982851E-12, 'S': -8.84655566777E-13, 't0':  '20041226.0060', 't1': '20060101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12060503696E-08, 'S':  5.40532443020E-09, 't0':  '20060101.0000', 't1': '20070101.0000'},
{'coef':'trnd', 'C': -2.23776998017E-12, 'S':  3.74692609570E-12, 't0':  '20060101.0000', 't1': '20070101.0000'},
{'coef':'acos', 'C':  4.26849040744E-14, 'S':  4.17013170120E-12, 't0':  '20060101.0000', 't1': '20070101.0000', 'per': 1.0},
{'coef':'asin', 'C': -3.11776943467E-12, 'S':  2.78580067808E-12, 't0':  '20060101.0000', 't1': '20070101.0000', 'per': 1.0},
{'coef':'acos', 'C': -8.14492693422E-13, 'S': -1.46128093369E-13, 't0':  '20060101.0000', 't1': '20070101.0000', 'per': 0.5},
{'coef':'asin', 'C':  1.34438187511E-12, 'S': -7.83684373020E-13, 't0':  '20060101.0000', 't1': '20070101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12082881396E-08, 'S':  5.40907135630E-09, 't0':  '20070101.0000', 't1': '20080101.0000'},
{'coef':'trnd', 'C':  1.34639840841E-12, 'S':  4.37108059157E-12, 't0':  '20070101.0000', 't1': '20080101.0000'},
{'coef':'acos', 'C':  5.65659322135E-13, 'S':  4.16060701415E-12, 't0':  '20070101.0000', 't1': '20080101.0000', 'per': 1.0},
{'coef':'asin', 'C':  8.93530262270E-13, 'S':  3.67148419721E-12, 't0':  '20070101.0000', 't1': '20080101.0000', 'per': 1.0},
{'coef':'acos', 'C': -1.31171281498E-12, 'S': -2.22007738616E-13, 't0':  '20070101.0000', 't1': '20080101.0000', 'per': 0.5},
{'coef':'asin', 'C':  2.17536978635E-12, 'S':  6.11044555895E-13, 't0':  '20070101.0000', 't1': '20080101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12069417412E-08, 'S':  5.41344243689E-09, 't0':  '20080101.0000', 't1': '20090101.0000'},
{'coef':'trnd', 'C': -4.58812450482E-12, 'S': -3.70704848481E-12, 't0':  '20080101.0000', 't1': '20090101.0000'},
{'coef':'acos', 'C': -6.14393797115E-13, 'S':  4.07771825537E-12, 't0':  '20080101.0000', 't1': '20090101.0000', 'per': 1.0},
{'coef':'asin', 'C': -9.58572764644E-13, 'S':  2.08522606461E-12, 't0':  '20080101.0000', 't1': '20090101.0000', 'per': 1.0},
{'coef':'acos', 'C': -7.23795469253E-14, 'S': -2.94226981352E-13, 't0':  '20080101.0000', 't1': '20090101.0000', 'per': 0.5},
{'coef':'asin', 'C':  6.11731219300E-13, 'S': -4.86992917127E-13, 't0':  '20080101.0000', 't1': '20090101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12115298657E-08, 'S':  5.40973538840E-09, 't0':  '20090101.0000', 't1': '20100227.0735'},
{'coef':'trnd', 'C':  4.14611092209E-12, 'S':  4.17696490922E-12, 't0':  '20090101.0000', 't1': '20100227.0735'},
{'coef':'acos', 'C': -4.42407095872E-14, 'S':  4.04236359433E-12, 't0':  '20090101.0000', 't1': '20100227.0735', 'per': 1.0},
{'coef':'asin', 'C': -8.25079543238E-13, 'S':  2.53126304987E-12, 't0':  '20090101.0000', 't1': '20100227.0735', 'per': 1.0},
{'coef':'acos', 'C': -7.16707788718E-13, 'S': -1.39866740064E-13, 't0':  '20090101.0000', 't1': '20100227.0735', 'per': 0.5},
{'coef':'asin', 'C':  7.95610838086E-13, 'S': -9.64905708157E-13, 't0':  '20090101.0000', 't1': '20100227.0735', 'per': 0.5},
{'coef':'gfct', 'C': -2.12099545763E-08, 'S':  5.41660118263E-09, 't0':  '20100227.0735', 't1': '20110311.0515'},
{'coef':'trnd', 'C': -3.15147196222E-12, 'S': -6.45495047703E-12, 't0':  '20100227.0735', 't1': '20110311.0515'},
{'coef':'acos', 'C':  1.76127614417E-12, 'S':  5.96470954145E-12, 't0':  '20100227.0735', 't1': '20110311.0515', 'per': 1.0},
{'coef':'asin', 'C':  2.26289275515E-13, 'S':  9.11159520381E-13, 't0':  '20100227.0735', 't1': '20110311.0515', 'per': 1.0},
{'coef':'acos', 'C': -1.77729899830E-12, 'S': -2.94457333992E-13, 't0':  '20100227.0735', 't1': '20110311.0515', 'per': 0.5},
{'coef':'asin', 'C':  1.22032432845E-12, 'S':  6.82898855183E-13, 't0':  '20100227.0735', 't1': '20110311.0515', 'per': 0.5},
{'coef':'gfct', 'C': -2.12082368477E-08, 'S':  5.41122033150E-09, 't0':  '20110311.0515', 't1': '20120101.0000'},
{'coef':'trnd', 'C': -6.11078657304E-12, 'S':  1.37916801459E-12, 't0':  '20110311.0515', 't1': '20120101.0000'},
{'coef':'acos', 'C':  1.85658015491E-12, 'S':  4.47092855258E-12, 't0':  '20110311.0515', 't1': '20120101.0000', 'per': 1.0},
{'coef':'asin', 'C': -3.15089651796E-12, 'S':  3.29241459498E-12, 't0':  '20110311.0515', 't1': '20120101.0000', 'per': 1.0},
{'coef':'acos', 'C': -1.38436765398E-12, 'S':  8.79315203250E-13, 't0':  '20110311.0515', 't1': '20120101.0000', 'per': 0.5},
{'coef':'asin', 'C':  2.31237125242E-12, 'S':  8.38302042916E-13, 't0':  '20110311.0515', 't1': '20120101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12131887847E-08, 'S':  5.41233795409E-09, 't0':  '20120101.0000', 't1': '20130101.0000'},
{'coef':'trnd', 'C':  5.21404107377E-12, 'S':  1.45220696672E-12, 't0':  '20120101.0000', 't1': '20130101.0000'},
{'coef':'acos', 'C': -1.60568790828E-13, 'S':  4.51663819090E-12, 't0':  '20120101.0000', 't1': '20130101.0000', 'per': 1.0},
{'coef':'asin', 'C':  7.89324573808E-14, 'S':  2.39640916963E-12, 't0':  '20120101.0000', 't1': '20130101.0000', 'per': 1.0},
{'coef':'acos', 'C':  6.00417674796E-13, 'S':  8.69188683256E-13, 't0':  '20120101.0000', 't1': '20130101.0000', 'per': 0.5},
{'coef':'asin', 'C':  1.78109739658E-12, 'S': -9.77867980041E-13, 't0':  '20120101.0000', 't1': '20130101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12079747436E-08, 'S':  5.41379016106E-09, 't0':  '20130101.0000', 't1': '20140101.0000'},
{'coef':'trnd', 'C': -7.86485104996E-12, 'S': -5.98437553842E-12, 't0':  '20130101.0000', 't1': '20140101.0000'},
{'coef':'acos', 'C':  9.20944257161E-13, 'S':  4.95000205108E-12, 't0':  '20130101.0000', 't1': '20140101.0000', 'per': 1.0},
{'coef':'asin', 'C': -7.58598253004E-12, 'S':  2.30009335917E-12, 't0':  '20130101.0000', 't1': '20140101.0000', 'per': 1.0},
{'coef':'acos', 'C': -4.50717355042E-13, 'S':  7.91241823090E-14, 't0':  '20130101.0000', 't1': '20140101.0000', 'per': 0.5},
{'coef':'asin', 'C': -5.20912301483E-13, 'S': -2.42716352743E-12, 't0':  '20130101.0000', 't1': '20140101.0000', 'per': 0.5},
{'coef':'gfct', 'C': -2.12158395947E-08, 'S':  5.40780578552E-09, 't0':  '20140101.0000', 't1': '20140615.0917'},
{'coef':'trnd', 'C':  2.46502884140E-12, 'S':  8.43773652778E-12, 't0':  '20140101.0000', 't1': '20140615.0917'},
{'coef':'acos', 'C':  1.68400415782E-12, 'S':  5.03233838519E-12, 't0':  '20140101.0000', 't1': '20140615.0917', 'per': 1.0},
{'coef':'asin', 'C':  3.18619810237E-12, 'S':  5.99671283004E-12, 't0':  '20140101.0000', 't1': '20140615.0917', 'per': 1.0},
{'coef':'acos', 'C': -1.25925099849E-12, 'S': -1.57870628503E-13, 't0':  '20140101.0000', 't1': '20140615.0917', 'per': 0.5},
{'coef':'asin', 'C':  1.12578755536E-12, 'S':  3.30587514133E-13, 't0':  '20140101.0000', 't1': '20140615.0917', 'per': 0.5},
{'coef':'gfct', 'C': -2.12147226548E-08, 'S':  5.41162904465E-09, 't0':  '20140615.0917', 't1': '20500101.0000'},
{'coef':'trnd', 'C':  0.00000000000E+00, 'S':  0.00000000000E+00, 't0':  '20140615.0917', 't1': '20500101.0000'},
{'coef':'acos', 'C':  3.04472200919E-13, 'S':  4.34034097939E-12, 't0':  '20140615.0917', 't1': '20500101.0000', 'per': 1.0},
{'coef':'asin', 'C': -1.53913671236E-12, 'S':  2.57940752999E-12, 't0':  '20140615.0917', 't1': '20500101.0000', 'per': 1.0},
{'coef':'acos', 'C': -1.06616774764E-12, 'S': -5.25926073450E-14, 't0':  '20140615.0917', 't1': '20500101.0000', 'per': 0.5},
{'coef':'asin', 'C':  1.33608393287E-12, 'S': -1.46418050564E-13, 't0':  '20140615.0917', 't1': '20500101.0000', 'per': 0.5}
]

def CS1608(t):
    C = 0e0
    S = 0e0
    for d in eigen_6s4:
        #print(d)
        t0 = d2s(d['t0'])
        t1 = d2s(d['t1'])
        #print('from {:} to {:}'.format(t0,t1))
        if t >= t0 and t < t1:
            if d['coef'] == 'gfct':
                C += d['C']
                #print('{:+.15e} {:.5e}'.format(d['C'], C))
                S += d['S']
            elif d['coef'] == 'trnd':
                diff = t - t0
                dt = (diff.days + diff.seconds/86400e0) / 365.25e0
                C += d['C'] * dt
                #print('{:+.15e} * {:.12e} {:.5e}'.format(d['C'], dt, C))
                S += d['S'] * dt
            elif d['coef'] == 'asin':
                diff = t - t0
                dt = (diff.days + diff.seconds/86400e0) / 365.25e0
                per = d['per']
                C += d['C'] * math.sin((math.pi * 2e0 / per)*dt)
                #print('{:+.15e} * {:.12e} {:.5e}'.format(d['C'], math.sin((math.pi * 2e0 / per)*dt), C))
                S += d['S'] * math.sin((math.pi * 2e0 / per)*dt)
            elif d['coef'] == 'acos':
                diff = t - t0
                dt = (diff.days + diff.seconds/86400e0) / 365.25e0
                per = d['per']
                C += d['C'] * math.cos((math.pi * 2e0 / per)*dt)
                #print('{:+.15e} * {:.12e} {:.5e}'.format(d['C'], math.cos((math.pi * 2e0 / per)*dt),C))
                S += d['S'] * math.cos((math.pi * 2e0 / per)*dt)
    return C, S

for year in range(1995, 2025):
    c,s = CS1608(datetime.datetime(year,1,1))
    print("C = {:+.15e} S = {:+.15e}".format(c,s))