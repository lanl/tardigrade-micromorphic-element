import numpy as np
import sympy

"""
tensor_tools.py

Tools for symbolic tensor manipulation

"""

def form_tensor(symbol,shape,isreal=True,isfunction=False,add_open_brace=False,add_close_brace=False,add_braces=False,variables=[]):
    
    if(add_braces):
        add_open_brace = True
        add_close_brace = True
    
    if(add_open_brace):
        symbol += r"{"
    
    if(len(shape)>1):
        return [form_tensor(symbol+str(i+1),shape[1:],\
                            isreal=isreal,isfunction=isfunction,\
                            add_close_brace=add_close_brace,variables=variables)\
                for i in range(shape[0])]
    else:
        
        if add_close_brace:
            cb = r"}"
        else:
            cb = ""
        
        if(isfunction):
            return [sympy.Function(symbol+str(i+1)+cb,real=isreal)(*variables) for i in range(shape[0])]
        else:
            return [sympy.symbols(symbol+str(i+1)+cb,real=isreal) for i in range(shape[0])]

def get_voigt_mapping(n):
    """
    Get the voigt index to true index mapping for a given number of indices.
    """
    
    voigt_indices = [[0,0],[1,1],[2,2],[1,2],[0,2],[0,1],[2,1],[2,0],[1,0]]
    
    if (n%2):
        indices = [[0],[1],[2]]
        n-=1
    else:
        indices = [[]]
        
    for _ in range(int(n/2)):
        
        temp = []
        
        for indx,i in enumerate(indices):
            
            for v in voigt_indices:
                temp.append(i+v)
                
        indices = temp
    return indices
    
def indices_to_matrix(in_string,tensor_name,index_sort):
    """
    Remove the indices from a string representing a tensor.
    
    index_sort = [number of indices in rows, number of indices in columns]
    
    order = sum(index_sort)
    
    """
    
    string = str(in_string)
    
    order = sum(index_sort)
    
    Ihats = get_voigt_mapping(index_sort[0])
    Jhats = get_voigt_mapping(index_sort[1])
    
    for i,Ihat in enumerate(Ihats):
        for j,Jhat in enumerate(Jhats):
            indices = ''
            a = []
            [a.append(v) for v in Ihat]
            [a.append(v) for v in Jhat]
            
            tensor_indices = ''.join([str(v+1) for v in  a])
            matrix_indices = r'(' + ','.join([str(i),str(j)]) + r')'
            
            string = string.replace(tensor_name + tensor_indices, tensor_name + matrix_indices)
    return string