
import numpy as np

def find_gap(a,b):
    position = []
    i = 0
    j = 0
    while(j < len(b)):
        if i >= len(a):
            position.append(i)
            j = j + 1
            continue

        if a[i] == b[j]:
            i = i + 1
            j = j + 1
        elif a[i] != b[j]:
            position.append(i)
            j = j + 1

    return position


def fill_gap(a,position):
    
    acc = 0
    for i in range(len(position)):
        a.insert(position[i]+acc,np.nan)
        acc += 1
    
    return a


