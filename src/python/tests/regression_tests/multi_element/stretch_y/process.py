import numpy as np
import os
import sys
sys.path.insert(0, '../../../..')
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt

import process_output

#Define extraction of value
def extract_values(UNSYMM,is_symm=False):

    #Form the unsymmetric tensor
    a = np.zeros([3,3])
    a[0,0] = UNSYMM[0]
    a[1,1] = UNSYMM[1]
    a[2,2] = UNSYMM[2]
    a[1,2] = UNSYMM[3]
    a[0,2] = UNSYMM[4]
    a[0,1] = UNSYMM[5]
    if(is_symm):
        a[2,1] = a[1,2]
        a[2,0] = a[0,2]
        a[1,0] = a[0,1]
    else:
        a[2,1] = UNSYMM[6]
        a[2,0] = UNSYMM[7]
        a[1,0] = UNSYMM[8]

    #Get the symmetric part
    a_symm   = 0.5*(a + a.T)
    a_unsymm = 0.5*(a - a.T)

    #Find the eigenvalues of a_symm

    w,v = np.linalg.eig(a_symm)

    #Get the axial vector
    a_vec = np.array([[-a_unsymm[1,2]],\
                      [ a_unsymm[0,2]],\
                      [-a_unsymm[0,1]]])

    #Rotate the axial vector to the principal basis
    return w,np.dot(v.T,a_vec)


dir = "../../../../../abaqus/tests/multi_element/stretch_y/"

response = process_output.ProcessMicroMorphic(dir)

disp   = np.linspace(0.,0.1,len(response.elements[1].steps[1].increments.keys())+1)

for el_num in response.elements.keys():
    PK2    = [[0]+[response.elements[el_num].steps[1].increments[v]['PK2'][w] for v in range(1,len(response.elements[1].steps[1].increments.keys())+1)] for w in range(9)]
    SIGMA  = [[0]+[response.elements[el_num].steps[1].increments[v]['SIGMA'][w] for v in range(1,len(response.elements[1].steps[1].increments.keys())+1)] for w in range(6)]
    labels = ['11','22','33','23','13','12','32','31','21']

    plt.figure()
    [plt.plot(disp,b,'-', label=c) for b,c in zip(PK2,labels)]
    [plt.plot(disp,b,'--',label=c) for b,c in zip(SIGMA,labels)]

    plt.xlabel(r"Deflection")
    plt.ylabel(r"PK2/SIGMA")
    plt.xlim([ 0.,0.15])
    plt.ylim([-0.2e9,1e9])
    plt.legend()
    plt.savefig('PK2_stretch_y_{0}.pdf'.format(el_num))

PC,AV = zip(*[extract_values(response.elements[1].steps[1].increments[v]['PK2'],is_symm=False) for v in range(1,len(response.elements[1].steps[1].increments.keys())+1)])

print PC[0],AV[0]
print PC[-1],AV[-1]
