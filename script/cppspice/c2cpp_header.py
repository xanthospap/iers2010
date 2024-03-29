#! /usr/bin/python

import sys, os

incdir = sys.argv[1]

head = """#ifdef __cplusplus
extern "C" {
#endif
"""

tail = """    #ifdef __cplusplus
}
#endif
"""

for file in os.listdir(incdir):
    if file.endswith(".h"):
        scratch = '.tmp.h'
        with open(scratch, 'w') as fout:
            print(head, file=fout)
            with open(os.path.join(incdir, file), 'r') as fin:
                for line in fin.readlines():
                    print(line.strip(), file=fout)
            print(tail, file=fout)
            os.rename(os.path.join(incdir, file), os.path.join(incdir, file+'.cvers'))
            os.rename(scratch, os.path.join(incdir, file))
