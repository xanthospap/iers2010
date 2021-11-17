#! /usr/bin/python
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

f1l=[]
f2l=[]

with open(f1) as fin:
  for line in fin.readlines():
    f1l.append(line)
with open(f2) as fin:
  for line in fin.readlines():
    f2l.append(line)

if len(f1l) != len(f2l):
  print("[ERROR] The two files do not have the same number of lines!", file=sys.stderr)
  sys.exit(1)

cols = ["mfh", "mfw"]

dif = {}
for c in cols: dif[c] = {'max': 0e0, 'sum': 0e0}

for rit, l in enumerate(f1l):
  l1 = l.split()
  l2 = f2l[rit].split()
  if len(l1) != len(cols):
    print('Line has {} cols but it should have {:}'.format(len(l1), len(cols)))
    sys.exit(2)
  if len(l2) != len(cols):
    print('Line has {} cols but it should have {:}'.format(len(l2), len(cols)))
    sys.exit(3)
  for cit, c in enumerate(cols):
    v1 = float(l1[cit])
    v2 = float(l2[cit])
    d = v1-v2
    dif[c]['sum'] += d
    if dif[c]['max'] < abs(d): 
      dif[c]['max'] = abs(d)

for k in dif:
  print('{:8s}: Max: {:.2e} Mean {:.2e}'.format(k, dif[k]['max'], dif[k]['sum']/len(f1l)))
