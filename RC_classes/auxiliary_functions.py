''' IMPORTING MODULES'''

import os
from warnings import warn

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions

'''
Some auxiliary functions for the tool
'''

def wrn(message):
    '''
    This function takes a string and print it in the log file
    
    Parameters
        ----------
        message : string
            message to be printed in the log file
    '''
    
    if not isinstance(message, str):
        try:
            message = str(message)
        except:
            raise TypeError(f'wrnn, message must be a string: message {message}')
    
    if not os.path.isdir(os.path.join('.','OutputReport')):
         os.mkdir('OutputReport')    
    output_path = os.path.join('.','OutputReport')
    
    if not os.path.isfile(os.path.join(output_path,'warnings.txt')):
         f = open(os.path.join(output_path,'warnings.txt'),"w")
         f.close()
         
    with open(os.path.join(output_path,'warnings.txt'),"a") as f:
        f.write(message + ' ')
    warn(message)