#! /usr/bin/python

import sys

hf = sys.argv[1]

head = """#ifdef __cplusplus
extern "C" {
#endif
"""

tail = """    #ifdef __cplusplus
}
#endif
"""

print(head)
with open(sys.argv[1], 'r') as fin:
    for line in fin.readlines():
        print(line.strip())
print(tail)
