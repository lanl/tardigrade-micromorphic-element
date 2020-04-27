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