import numpy as np
import Bio.PDB as bpdb

"""Subprogram for generating water orientations"""

def generate_new_water_orientation(w, threshold):
    "This generates a new water orientation starting from a water residue with RSCC scores calculated"
    d1 = w['D1'].get_vector()
    d2 = w['D2'].get_vector()
    o = w['O'].get_vector()
    rotd1 = d1
    rotd2 = d2
     
    d1 = d1 - o
    d2 = d2 - o
    

    if(w['D1'].bfactor >= (threshold - 0.2) and w['D1'].bfactor < threshold):
        step1 = 5
        m1 = bpdb.vectors.rotaxis(np.pi*step1/180, d2)
        rotd1 = d1.left_multiply(m1)
        w['D1'].set_coord(rotd1._ar+o._ar)
    elif(w['D1'].bfactor < (threshold - 0.2)):
        step1 = 10
        m1 = bpdb.vectors.rotaxis(np.pi*step1/180, d2)
        rotd1 = d1.left_multiply(m1)
        w['D1'].set_coord(rotd1._ar+o._ar)

        
    if(w['D2'].bfactor >= (threshold - 0.2) and w['D2'].bfactor < threshold):
        step2 = 5
        m2 = bpdb.vectors.rotaxis(np.pi*step2/180, rotd1)
        rotd2 = d2.left_multiply(m2)
        w['D2'].set_coord(rotd2._ar+o._ar)
    elif(w['D2'].bfactor < (threshold - 0.2)):
        step2 = 10
        m2 = bpdb.vectors.rotaxis(np.pi*step2/180, rotd1)
        rotd2 = d2.left_multiply(m2)
        w['D2'].set_coord(rotd2._ar+o._ar)
        
    return w

