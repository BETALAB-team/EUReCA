'''IMPORTING MODULES'''

import pandas as pd
import numpy as np
from scipy import interpolate
from RC_classes.auxiliary_functions import wrn

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions

def loadEnvelopes(path):
    
    '''
    loadEnvelopes loads all materials, windows, stratigraphies and envelopes
    archetypes in a database that is then utilized by the main program
    see file \input\Buildings_Envelopes.xlsx
    
    Parameters
    ----------
    path : string
        Path containing the string of the Buildings_Envelopes.xlsx

    Returns
    -------
    envelopesDict: dictionary with archetype_key/Envelope(object) data
    
    '''
    
    # Check input data type
    
    if not isinstance(path, str):
        raise TypeError(f'Ops... input path is not a string: path {path}') 
    
    # Tries to import the 4 excel sheet with: meterial, windows, stratigraphies, envelopes
    # Creates 4 dataframes with this tables
    
    try:
        materials = pd.read_excel(path,sheet_name="Materials",header=[1],index_col=[0])
        windows = pd.read_excel(path,sheet_name="Windows",header=[1],index_col=[0])
        stratigraphies = pd.read_excel(path,sheet_name="Constructions",header=[0],index_col=[0])
        envelopes = pd.read_excel(path,sheet_name="Envelopes",header=[0],index_col=[0])
    except FileNotFoundError:
        raise FileNotFoundError(f'Failed to open the archetype xlsx file {path}... Insert a proper path')
    
    # MATERIALS:    the list of materials is cycled and for each of them a material object is created
    #               the materials objects are stored in a dictionary materialsDict
    
    materials = materials[1:]
    materialsDict=dict()
    for i in materials.index:
        material = OpaqueMaterial(materials.loc[i])
        materialsDict[i]=material
    
    # WINDOWS:      the list of windows is cycled and for each of them a windows object is created
    #               the windows objects are stored in a dictionary windowsDict
       
    windowsDict=dict()
    for i in windows.index:
        window = Window(windows.loc[i])
        windowsDict[i]=window
    
    # STRATIGRAPHIES:   the list of Structures is cycled and for each of them a OpaqueStratigraphy object is created
    #                   the Structures objects are stored in a dictionary stratigraphiesDict
    #                   In order to be created, the OpaqueStratigraphy object needs both the list of materials it is made, and the stratigraphiesDict
    
    stratigraphiesDict=dict()
    for i in stratigraphies.index:
        stratigraphy = OpaqueStratigraphy(stratigraphies.loc[i],materialsDict)
        stratigraphiesDict[i]=stratigraphy
    
    # ENVELOPES:    the list of Envelopes is cycled and for each of them a Envelope object is created
    #               the Envelopes objects are stored in a dictionary envelopesDict
    #               In order to be created, the Envelope object needs both the list of Structures it is made, both the stratigraphiesDict, and the windowsDict
    
    envelopesDict=dict()
    for i in envelopes.index:
        envelope = Envelope(envelopes.loc[i],stratigraphiesDict,windowsDict)
        envelopesDict[envelope.name]=envelope

    return envelopesDict

#%%--------------------------------------------------------------------------------------------------- 
#%% Envelope class

class Envelope(object):
    '''
    Definition of the Envelope class
    Each object of Envelope type contains several informations about stratigraphies
    
    __init__ method takes:
        envelopeslist: dictionary containing all stratigraphies codes
        stratigraphiesDict: dictionary with all the stratigraphies
        windowsDict: dictionary with all the windows
        
    Methods:
        init
    '''

    def __init__(self,envelopeList,stratigraphiesDict,windowsDict):
        '''
        Initializes an envelope object
        
        Parameters
            ----------
            envelopeList : pandas series
                pandas series describing the list of keys of each stratigraphy
            stratigraphiesDict : dictionary
                dictionary with all the stratigraphies objects with their code as key
            windowsDict : dictionary
                dictionary with all the windows objects with their code as key
                
        Returns
        -------
        None.
        
        '''  

        # Check input data type
        
        if not isinstance(envelopeList, pd.core.series.Series):
            raise TypeError(f'Ops... input envelopeList is not a pandas series: envelopeList {envelopeList}') 
        if not isinstance(stratigraphiesDict, dict):
            raise TypeError(f'Ops... input stratigraphiesDict is not a dictionary: stratigraphiesDict {stratigraphiesDict}') 
        if not isinstance(windowsDict, dict):
            raise TypeError(f'Ops... input windowsDict is not a dictionary: windowsDict {windowsDict}') 
    
        try:
            self.id = envelopeList.index
            self.name = envelopeList["name"]
            self.Roof = stratigraphiesDict[envelopeList["Roof"]]
            self.GroundFloor = stratigraphiesDict[envelopeList["GroundFloor"]]
            self.IntCeiling = stratigraphiesDict[envelopeList["IntCeiling"]]
            self.ExtWall = stratigraphiesDict[envelopeList["ExtWall"]]
            self.IntWall = stratigraphiesDict[envelopeList["IntWall"]]
            self.IntFloor = stratigraphiesDict[np.flip(envelopeList["IntCeiling"])]                                                                    
            self.Window = windowsDict[envelopeList["Window"]]
        except KeyError:
            raise KeyError(f"There's something wrong with the Envelopes creation, check the envelopes sheet in Buildings_Envelopes.xlsx")

#%%--------------------------------------------------------------------------------------------------- 
#%% OpaqueMaterial class

class OpaqueMaterial(object):
    
    '''
    Definition of the class opaqueMaterial
    Contains all materials properties
    
    __init__ method takes:
        PropList: list of the properties of the material
        
    Methods:
        init
    '''

    def __init__(self,PropList):
        '''
        Initializes an OpaqueMaterial object
        
        Parameters
            ----------
            PropList: pandas series
                pandas series describing the list of materials properties
            
        Returns
        -------
        None.
        
        '''  
        # Check input data type
        
        if not isinstance(PropList, pd.core.series.Series):
            raise TypeError(f'Ops... input PropList is not a pandas series: PropList {PropList}') 
               
        try:
            self.name = PropList['name']
            self.type = PropList['Material']
            self.s = PropList['thickness']
            self.l = PropList['conductivity']
            self.d = PropList['density']
            self.c = PropList['specific_heat']
            self.a = PropList['thermal absorptance']                               # Shortwave for vdi-6007
        except KeyError:
            raise KeyError(f"There's something wrong with the Materials creation, check the materials sheet in Buildings_Envelopes.xlsx")
        
        # Check materials properties
        
        if self.s < 0.001 or self.s > 1.:
            wrn(f"\n\nMaterial {self.name}. Are you sure about the material thickness?? thickness {self.s} m\n")
        if self.l < 0.01 or self.l > 10.:
            wrn(f"\n\nMaterial {self.name}. Are you sure about the material conductivity?? conductivity {self.l} W/(m K)\n")
        if self.d < 1. or self.d > 10000.:
            wrn(f"\n\nMaterial {self.name}. Are you sure about the material density?? density {self.d} kg/m3")
        if self.c < 500. or self.c > 2000.:
            wrn(f"\n\nMaterial {self.name}. Are you sure about the material specific heat?? specific heat {self.c} J/(kg K)\n")
        if self.a < 0. or self.a > 1.:
            wrn(f"\n\nMaterial {self.name}. Are you sure about the material absorption coefficient?? absorption coefficient {self.a}\nThe coefficient is set to 0.5\n")
            self.a = 0.5
               
        self.massless = 0

        if self.type == "NoMass":
            self.r = PropList['Resistance']
            self.massless = self.r
        else:
            self.r = self.s/self.l
            
        if self.r > 5.:
            wrn(f"\n\nMaterial {self.name}. The thermal resistance is pretty high.. thermal resistance {self.r} (m2 K)/W\n")
                   

#%%--------------------------------------------------------------------------------------------------- 
#%% OpaqueMaterial class

class OpaqueStratigraphy(object):
    
    '''
    Defines stratigraphies class
    and calculates all the stratigraphy parameter
    U,k_int,k_est
    
    __init__:
        stratigraphy: list of materials codes
        materialsDict: dictionary of all the materials
        
    ISO13790params calculates the params for the 1C thermal network: no input requested
    
    vdi6007params calculates the params for the 2C thermal network: no input requested
    
    vdi6007surfaceParams calculates some params for a specific surface:
        area: area of the surface [m2]
        Asim: flag to indicate if the surface is either assimmetric or not (vdi). boolean True or False
    
    Methods:
        init
        ISO13790params
        vdi6007params
        vdi6007surfaceParams
        printInfo
    '''
    
    # Class attributes
    
    alpha_lim = pd.DataFrame({'Outside':[25,25,1000,7.7,7.7],
                             'Inside':[7.7,7.7,7.7,7.7,7.7]},
                            index=['ExtWall','Roof','GroundFloor','IntWall','IntCeiling'])
    
    alpha_rad = 5
    
    def __init__(self,stratigraphy,materialsDict):
        '''
        Initializes an OpaqueStratigraphy object
        
        Parameters
            ----------
            stratigraphy: pandas series
                pandas series describing the list of materials the structure is made off
            materialsDict: dictionary
                dictionary including all materials objects
            
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        if not isinstance(stratigraphy, pd.core.series.Series):
            raise TypeError(f'Ops... input stratigraphy is not a pandas series: stratigraphy {stratigraphy}') 
        if not isinstance(materialsDict, dict):
            raise TypeError(f'Ops... input materialsDict is not a dictionary: materialsDict {materialsDict}') 
               
        try:        
            self.name = stratigraphy['name']
            self.wallType = stratigraphy['opaque_type_name']
            self.listmaterials = stratigraphy
            self.listmaterials = [x for x in stratigraphy if str(x) != 'nan']
            self.listmaterials = self.listmaterials[4:]
        except KeyError:
            raise KeyError(f"There's something wrong with the Stratigraphies creation, check the stratigraphies sheet in Buildings_Envelopes.xlsx")
        
        # Set outside and inside convective and radiant heat transfer coefficients
                
        self.R_si = 1/(self.alpha_lim.loc[self.wallType]['Inside'])
        self.R_se = 1/(self.alpha_lim.loc[self.wallType]['Outside'])
        self.alpha_conv_int = self.alpha_lim.loc[self.wallType]['Inside'] - self.alpha_rad
        self.alpha_conv_est = self.alpha_lim.loc[self.wallType]['Outside'] - self.alpha_rad
        
        # Pick ups materials from the dictionary
        
        try:
            for i in range(len(self.listmaterials)):
                self.listmaterials[i] = materialsDict[self.listmaterials[i]]              
        except KeyError:
            raise KeyError(f"There's something wrong with the {self.name} Stratigraphy creation, I can not find all materials.  Check the stratigraphies and materials sheets in Buildings_Envelopes.xlsx")
        
        # Creates a matrix with the properties of the layers
        
        self.alpha_est = self.listmaterials[0].a                               # From the external to the internal side
        self.nL=len(self.listmaterials)
        self.properties=np.zeros([self.nL,8])
        for i in range(self.nL):
            self.properties[i,0]=self.listmaterials[i].s
            self.properties[i,1]=self.listmaterials[i].l
            self.properties[i,2]=self.listmaterials[i].d
            self.properties[i,3]=self.listmaterials[i].c
            self.properties[i,4]=self.listmaterials[i].r
            self.properties[i,5]=self.listmaterials[i].massless
        
        # Calculation of the U-valus [W/(m2 K)]
        
        self.U_net = 1/(sum(self.properties[:,4]))
        self.U = 1/(sum(self.properties[:,4])+self.R_si+self.R_se)
        
        if self.wallType == 'GroundFloor':
            self.U=self.U*0.7
            
        # Parameters of the structure
        
        self.s = np.array(self.properties[:,0])
        self.cond = np.array(self.properties[:,1])
        self.rho = np.array(self.properties[:,2])
        self.cp = np.array(self.properties[:,3])
        self.r = np.array(self.properties[:,4])
        self.massless = np.array(self.properties[:,5])
        
        # Run ISO13790params and vdi6007params to calculate further parameters
        
        self.ISO13790params()
        self.vdi6007params()


    def ISO13790params(self):
        '''
        Calculates ISO13790 params
        
        Parameters
            ----------
            None.
            
        Returns
        -------
        None.
        '''    

        # Set some parameters from the standard
        
        T=86400
        sigma_2 = (T/np.pi)*(self.cond/(self.rho*self.cp))
        sigma = np.sqrt(sigma_2)                                               # Depth of penetration
        eps = self.s/sigma

        # Thermal transfer matrix
        
        Z = np.zeros((2,2,self.nL),complex)

        for i in range(self.nL):

            Z[0,0,i] = np.cosh(eps[i])*np.cos(eps[i])+1j*np.sinh(eps[i])*np.sin(eps[i])
            Z[1,1,i] = Z[0,0,i]
            Z[0,1,i] = -(sigma[i]/(2*self.cond[i]))*(np.sinh(eps[i])*np.cos(eps[i])+np.cosh(eps[i])*np.sin(eps[i]) + 1j*(np.cosh(eps[i])*np.sin(eps[i])-np.sinh(eps[i])*np.cos(eps[i])))
            Z[1,0,i] = -(self.cond[i]/sigma[i])*(np.sinh(eps[i])*np.cos(eps[i])-np.cosh(eps[i])*np.sin(eps[i]) + 1j*(np.sinh(eps[i])*np.cos(eps[i])+np.cosh(eps[i])*np.sin(eps[i])))

        Z_si = np.eye(2)
        Z_si[0,1] = -self.R_si                                                 # Internal surface resistance (convection and radiation, ISO 6946)
        Z_se = np.eye(2)
        Z_se[0,1] = -self.R_se                                                 # External surface resistance (convection and radiation, ISO 6946)
        i=self.nL-1

        while i > 0:
            Z[:,:,i-1] = np.matmul(Z[:,:,i],Z[:,:,i-1])
            i=i-1

        Z = Z[:,:,0]
        Z = np.matmul(Z_se,Z)
        Z = np.matmul(Z,Z_si)                                                  # Thermal transfer matrix

        # Calculation of the dynamic thermal characteristics
        
        R_w_stat = self.R_se + self.R_si + np.sum(self.r)                      # Global thermal resistance [(m2 K)/W]
        U_0 = 1/R_w_stat

        # Thermal admittances (2) and specific thermal capacity (2)
        
        #Y = np.zeros((2,2),complex)

        #Y[0,0] = -Z[0,0]/Z[0,1]
        #Y[1,1] = -Z[1,1]/Z[0,1]
        #Y_11 = np.abs(Y[0,0])                                                 # (1,1) internal side
        #Y_22 = np.abs(Y[1,1])                                                 # (2,2) external side
        #deTy_11 = T/(2*np.pi)*np.arctan2(np.imag(Y[0,0]),np.real(Y[0,0]))/3600  # Temporal variation of admittances [h]
        #deTy_22 = T/(2*np.pi)*np.arctan2(np.imag(Y[1,1]),np.real(Y[1,1]))/3600  # Temporal variation admittances
        self.k_est = (T/(2*np.pi))*np.abs((Z[0,0]-1)/Z[0,1])                   # Thermal capacitance
        self.k_int = (T/(2*np.pi))*np.abs((Z[1,1]-1)/Z[0,1])                   # Thermal capacitance
        self.k_mean = (self.k_int + self.k_est)/2
        #Y_12 = np.abs(-1/Z[0,1])                                              # Trasmittanza termica periodica
        #f = Y_12/U_0                                                          # Fattore di decremento
        
                                                    
    def vdi6007params(self):
        '''
        Calculates vdi6007 params
        
        Section 6.3
        
        # vdi6007params Calculates the parameters (thermal resistance and thermal
        # capacitance) associated with the building envelope LP according to the 2-c
        # model of standard VDI 6007
        
        # Inputs
        #   Matrix LP describes building envelope; each row is a building element.
        #   Columns are thickness (s), thermal conductivity (cond), density (rho) 
        #   and specific heat (cp)
        #   Total surface area S of wall with building envelope LP
        #   Flag 'asim' indicates whether building component LP is asimmetrically
        #   loaded (asim = 1) or not (asim = 0) because in the first case C1_korr 
        #   must be considered instead of C1
        
        # Subscripts
        # AW external walls and internal walls facing unheated areas
        # IW internal walls
        
        # Period T_bt 
        # T_bt = 7 days for a single building component; 
        # T_bt = 2 days for building components where thermal storage masses are 
        # thermally covered on the room side (eg suspended ceilings)
        # Calculation are conducted with both periods, then the resulting
        # parameters R1 and C1 are compared and the right T_bt is chosen according
        # to the criterion given in the standard
        
        # Outputs
        #    R1 - dynamic rhermal resistance [K/W]
        #    C1 - dynamic thermal capacity [K/W]
        #    Rw - static specific thermal resistance [(m2 K) / W]
        
        # Determine thermal resistance R and thermal capacitance of layers
        
        Parameters
            ----------
            None.
            
        Returns
        -------
        None.
        '''
    
        R = self.r         # layers thermal resistance [m2 K / W]
        C = self.cp*self.rho*self.s      # layers thermal capacitance [J m2 / K]
        self.C = sum(C)               
    
        T_bt = np.array([2,7])    # days
        self.omega_bt = 2*np.pi/(86400*T_bt)
        #T_ra = 5        # days
        #omega_ra = 2*pi./(86400*T_ra) 
        
        Z_t2 = np.zeros((2,2,self.nL),complex) 
        Z_t7 = np.zeros((2,2,self.nL),complex) 
    
        for i in range(self.nL):
            
            r = R[i]
            c = C[i] 
            '''
            R1 = np.zeros((2,1))
            R2 = R1
            C1 = R1
            C2 = R1 
            '''
            #r_ah = self.massless[i]
            Av = np.zeros((2,2,2),complex)             
            
            #if r_ah == 0:
                      
            for om in range(2):
                arg = np.sqrt(0.5*self.omega_bt[om]*r*c)           

                Re_a11 = np.cosh(arg)*np.cos(arg)
                Im_a11 = np.sinh(arg)*np.sin(arg) 
                Re_a22 = Re_a11
                Im_a22 = Im_a11 
    
                Re_a12 = r/(2*arg)*           (   np.cosh(arg)*np.sin(arg)       +      np.sinh(arg)*np.cos(arg)       ) 
                Im_a12 = r/(2*arg)*           (   np.cosh(arg)*np.sin(arg)       -      np.sinh(arg)*np.cos(arg)       )
 


                Re_a21 = -arg/r*(np.cosh(arg)*np.sin(arg) - np.sinh(arg)*np.cos(arg)) 
                Im_a21 = arg/r*(np.cosh(arg)*np.sin(arg) + np.sinh(arg)*np.cos(arg)) 

                Av[0,0,om] = Re_a11 + 1j*Im_a11 
                Av[1,1,om] = Re_a22 + 1j*Im_a22 
                Av[0,1,om] = Re_a12 + 1j*Im_a12 
                Av[1,0,om] = Re_a21 + 1j*Im_a21 

            Z_t2[:,:,i] = Av[:,:,0] 
            Z_t7[:,:,i] = Av[:,:,1] 
                
            #else:
                
            #    Z_t2[:,:,i] = np.eye(2) 
            #    Z_t2[0,1,i] = -r_ah 
            #    Z_t7[:,:,i] = Z_t2[:,:,i]
                
            #    #R(i) = r_ah 
                
        # Chain matrix of the complete wall A_1n
        
        # for i in range(self.nL-1):
        #     Z_t2[:,:,i+1]  = np.matmul(Z_t2[:,:,i],Z_t2[:,:,i+1]) 
        #     Z_t7[:,:,i+1]  = np.matmul(Z_t7[:,:,i],Z_t7[:,:,i+1])
        # self.A1n_t2 = np.zeros((2,2,1),complex)
        # self.A1n_t7 = np.zeros((2,2,1),complex)
        # self.A1n_t2 = Z_t2[:,:,self.nL-1] 
        # self.A1n_t7 = Z_t7[:,:,self.nL-1]
        self.A1n_t2 = np.zeros((2,2,1),complex)
        self.A1n_t7 = np.zeros((2,2,1),complex)
        self.A1n_t2 = Z_t2[:,:,-1]
        self.A1n_t7 = Z_t7[:,:,-1]
        for t in range(-2,-self.nL-1,-1):
            self.A1n_t2  = np.matmul(self.A1n_t2,Z_t2[:,:,t])
            self.A1n_t7  = np.matmul(self.A1n_t7,Z_t7[:,:,t])

        # self.A1n_t2 = Z_t2[:,:,0] 
        # self.A1n_t7 = Z_t7[:,:,0]               
        
    
    def vdi6007surfaceParams(self,sup,asim):
        '''
        Calculates vdi6007 params
        
        Section 6.4
        
        Parameters
            ----------
            sup : float
                area of the surface.
            asim: boolean
                Is the surface non-adiabatic? True/False
                
        Returns
        -------
        R1,C1: tuple of floats
            Resistance and capacitance R1 and C1
        '''
        
        # Check input data type
        
        if not isinstance(sup, float):
            raise TypeError(f'Ops... surface {self.name} input sup is not a float: sup {sup}') 
        if not isinstance(asim, bool):
            raise TypeError(f'Ops... surface {self.name} input asim is not a boolean: asim {asim}') 

        # Procudeure Section 6.4
        
        if sup == 0:
            sup = 0.0000001
        rw = sum(self.r)/sup
        
        R1_t = dict()
        C1_t = dict()
        
        for a, omega, days in zip([self.A1n_t2,self.A1n_t7],[self.omega_bt[0],self.omega_bt[1]],['2','7']):
            #rcValues Given the complex matrix of the building element BT, the function
            #calculates the values R1 and C1        
            R1 = 1/sup*((np.real(a[1,1])-1)*np.real(a[0,1])+np.imag(a[1,1])* np.imag(a[0,1])  ) / ( ( np.real(a[1,1])-1 )**2 + (np.imag(a[1,1]))**2 )
            
            #R2 = 1/sup *(  (np.real(a [0,0])-1) * np.real(a [0,1]) + np.imag(a [0,0])*np.imag(a [0,1]) ) / ((np.real(a [0,0])-1)^2+np.imag(a [0,0])^2)
            
            if asim == False:
                C1 = sup*  ((np.real(a [1,1]) - 1)**2  +  (np.imag(a [1,1]))**2 )  /  (omega * (np.real(a [0,1])*np.imag(a [1,1]) - (np.real(a [1,1])-1)*np.imag(a [0,1])  )  )
            else:
                # sarebbe C1_korr per pareti caricate asimmetricamente (pareti AW)
                C1 = sup*( 1 / (omega*R1*sup)) * (rw*sup - np.real(a [0,1])*np.real(a [1,1]) - np.imag(a [1,1]) * np.imag(a [0,1]) ) / (np.real(a [1,1])*np.imag(a [0,1]) - np.real(a [0,1])*np.imag(a [1,1]) )

            R1_t[days] = R1
            C1_t[days] = C1
        
        rr = R1_t['2']/R1_t['7']
        cr = C1_t['2']/C1_t['7']
        
        # versione di jacopo
        #if (rr>0.99 and cr<0.95) or (((rr<0.99 and cr<0.95) and (abs(rr-cr)>0.3))):
        if (rr>0.99 and cr<0.95) or (((rr<0.95 and cr<0.95) and (abs(rr-cr)>0.3))):
            R1 = R1_t['2'] # T_bt = 2 days
            C1 = C1_t['2']
        else:
            R1 = R1_t['7'] # T_bt = 7 days
            C1 = C1_t['7']
            
        return R1, C1
    
    def printInfo(self):
        '''
        Just print some info about the structure
        
        Parameters
            ----------
            None
                
        Returns
        -------
        None
        '''
        print(self.name)
        print(self.U)
        print(self.k_int)
        print(self.k_est)

#%%--------------------------------------------------------------------------------------------------- 
#%% Window class

class Window(object):
    
    '''
    Defines the simple window model with all its characteristics
    and apply the simple glazing model to estimates the curve of SHGC(theta)
    
    __init__: 
        propList: list of the properties of the window
        
    simpleGlazingModel generates the profiles of SHGC. No input
    
    Methods:
        init
        simpleGlazingModel
    ''' 
    
    def __init__(self,propList):
        '''
        initialize Window object
        
        Parameters
        ----------
        propList : pandas series
                pandas series describing the list of the window properties    
        
        Returns
        -------
        None
        '''
        
        # Check input data type
        
        if not isinstance(propList, pd.core.series.Series):
            raise TypeError(f'Ops... input PropList is not a pandas series: PropList {PropList}') 
                
        try:
            self.id = propList.index
            self.name = str(propList["name"])
            self.U = float(propList["U"])
            self.SHGC = float(propList["SHGC"])
            self.Tv = float(propList["tau_vis"])
            self.F_f = float(propList["F_f"])
            self.F_sh =float(propList["F_sh"])
            self.F_so = float(propList["F_so"])
            self.F_w = float(propList["F_w"])  
            self.wwr = [float(propList["wwr_N"]),float(propList["wwr_E"]),float(propList["wwr_S"]),float(propList["wwr_W"])]
        except KeyError:
            raise KeyError(f"There's something wrong with the Window creation, check the windows sheet in Buildings_Envelopes.xlsx")
        
        # Check properties quality
        
        if self.U < 1. or self.U > 7.:
            wrn(f"\n\nWindow {self.name}. Are you sure about the window's U-value?? U-value {self.U} W/(m2 K)\n")
        if self.SHGC < 0.01 or self.SHGC > 1.:
            wrn(f"\n\nWindow {self.name}. Are you sure about the window's SHGC?? SHGC {self.SHGC}, SHGC is set to 0.7\n")
            self.SHGC = 0.7
        if self.Tv < 0.01 or self.Tv > 1.:
            wrn(f"\n\nWindow {self.name}. Are you sure about the window's visible transmittance?? Tv {self.Tv}, Tv is set to 0.7\n")
            self.Tv = 0.7
        if self.F_f < 0.0 or self.F_f > 1.:
            wrn(f"\n\nWindow {self.name}. Are you sure about the window's frame fraction?? F_f {self.F_f}, F_f is set to 0.7\n")
            self.F_f = 0.7
        if self.F_sh < 0.0 or self.F_sh > 1.:
            wrn(f"\n\nWindow {self.name}. Are you sure about the window's shading coefficient?? F_sh {self.F_sh}, F_sh is set to 0.7\n")
            self.F_sh = 0.7
        if self.F_so < 0.0 or self.F_so > 1.:
            wrn(f"\n\nWindow {self.name}. Are you sure about the window's F_so?? F_so {self.F_so}, F_so is set to 0.7\n")
            self.F_so = 0.7
        for w_w_r in self.wwr:
            if w_w_r > 1.:
                wrn(f"\n\nWindow {self.name}. Are you sure about the window's window to wall ratio?? w_w_r {self.w_w_r}\n")
            
        # Runs the simpleGlazingModel    
            
        self.simpleGlazingModel()
        
    def simpleGlazingModel(self):
        
        '''
        Description
        This function implement the simplified model for glazing that is used in energy plus
        The output variables are:
        -   4 vectors containing solar transmittance, reflectance, absorptance, SHGC 
            in function of the incident angle  from 0� to 90�, with a 10� step 
        -   the solar transmittance, reflectance and absorption at normal incident;
        -   the equivalent conductivity [W/m K] of the glazing material;
        -   the equivalent thickness [m] of the glazing layer;
        -   the visible reflectance on the front side and back side;
        -   the inward flowing fraction;
        
        ---------ATTENTION-----------
        For some combinations of variable SHGC and U the curves are not
        physically possible,(reflectance > 1). 
        For this reason it is necessary to manipulate the final reflectance and transmittance 
        vectors to manage the curves.
        It is also necessary to force the absolute reflectance to be 1 and the
        absolute transmittance to be 0 for the perpendicular incidence.
        
        
        Controlling input variable
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None        
        '''
        
        if self.U < 0:
            print('Negative window U-value. Simulation won''t proceed.')
            return
        
        if self.U > 7:
            print('U-value may be to high, the model could evaluate un-proper output variable.')
        
        if self.SHGC < 0 or self.SHGC >1:
            print('Solar gain of the window not allowed. Simulation won''t proceed.');
            return
        
        # Calculation of thermal resistances W/m^2 K
        
        if self.U < 5.85:                                
            self.Ri_w = 1/(0.359073*np.log(self.U)+6.949915)
            
        else:
            self.Ri_w = 1/( 1.788041*self.U - 2.886625)       
        
        self.Ro_w = 1/(0.025342*self.U + 29.163853)
        self.Rl_w = 1/self.U- self.Ro_w - self.Ri_w

        if 1/self.Rl_w > 7:                                                    # Thickness d of the equivalent layer, m 
            self.d = 0.002             
        else:
            self.d = 0.05914 - 0.00714/self.Rl_w             
        
        self.Keff = self.d/self.Rl_w                                           # Thermal conductivity of the representative layer, W/m K

        # Calculation of solar transmittance for normal incident
        
        if  self.U > 4.5:
            if self.SHGC < 0.7206:
                self.Ts = (0.939998 * self.SHGC**2) + (0.20332 * self.SHGC)    
            else:
                self.Ts= (1.30415*self.SHGC) - 0.30515
            
        elif self.U < 3.4:            
            if self.SHGC <= 0.15:
                self.Ts = 0.41040*self.SHGC 
            else:
                self.Ts = (0.085775*self.SHGC**2) + (0.963954*self.SHGC) - 0.084958     
            
        else:
            if self.SHGC <= 0.15:
                self.Ts_1 = (0.939998 * self.SHGC**2) + (0.20332 * self.SHGC) 
                self.Ts_2 = 0.41040*self.SHGC
                self.Ts = self.Ts_2 + (self.Ts_1-self.Ts_2)/(4.5-3.4)*(self.U-3.4)
            
            elif self.SHGC > 0.15 and self.SHGC < 0.7206:
                self.Ts_1 = (0.939998 * self.SHGC**2) + (0.20332 * self.SHGC)
                self.Ts_2 = (0.085775*self.SHGC**2) + (0.963954* self.SHGC) - 0.084958              
                self.Ts = self.Ts_2 + (self.Ts_1-self.Ts_2)/(4.5-3.4)*(self.U-3.4)
            
            else:
                self.Ts_1 = (1.30415*self.SHGC) - 0.30515
                self.Ts_2 = (0.085775*self.SHGC**2) + (0.963954* self.SHGC) - 0.084958
                self.Ts = self.Ts_2 + (self.Ts_1-self.Ts_2)/(4.5-3.4)*(self.U-3.4)
         
                
        # Calculation of inside and outside film resistances in summer condition 
        
        self.x = self.SHGC - self.Ts
        
        if self.U > 4.5:
            self.Ri_s = 1/((29.436546*self.x**3) - (21.943415*self.x**2) + (9.945872*self.x) + 7.426151)
            self.Ro_s = 1/((2.225824*self.x) + 20.577080)

        elif self.U < 3.4:
            self.Ri_s = 1/((199.8208128*self.x**3) - (90.639733*self.x**2) + (19.737055*self.x) + 6.766575)
            self.Ro_s = 1/((4.475553*self.x) + 20.674424)
        
        else:
            self.Ri_s_1 = 1/((29.436546*self.x**3) - (21.943415*self.x**2) + (9.945872*self.x) + 7.426151)
            self.Ri_s_2 = 1/((199.8208128*self.x**3) - (90.639733*self.x**2) + (19.737055*self.x) + 6.766575)
            self.Ri_s = self.Ri_s_2 + (self.Ri_s_1-self.Ri_s_2)/(4.5 - 3.4)*(self.U - 3.4)
            
            self.Ro_s_1 = 1/((2.225824*self.x) + 20.577080)
            self.Ro_s_2 = 1/((4.475553*self.x) + 20.674424)
            self.Ro_s = self.Ro_s_2 + (self.Ro_s_1-self.Ro_s_2)/(4.5 - 3.4)*(self.U - 3.4)
            
        # Inward flowing fraction
        self.Rl = self.Rl_w
        self.N = (self.Ro_s + 0.5*self.Rl)/(self.Ro_s+self.Rl+self.Ri_s)
        
        self.As = self.x/self.N
        self.Rs_f= 1 - self.Ts - self.As
        self.Rs_b = self.Rs_f
        
        # Evaluating the visible proprties of the equivalent layer
        
        self.Rv_f = -0.0622*self.Tv**3 + 0.4277*self.Tv**2 - 0.4169*self.Tv + 0.2399
        self.Rv_b = -0.7409*self.Tv**3 + 1.6531*self.Tv**2 - 1.2299*self.Tv + 0.4545
        
        # Definition of the polinomials 
        # Loading polynomial coefficients
        
        TransCoef = [[0.014700000,1.486000000,-3.852000000,3.355000000,-0.001474000],\
                     [ 0.554600000,	0.035630000,-2.416000000,2.831000000,-0.002037000],
                     [0.770900000,	-0.638300000,	-1.576000000,	2.448000000,	-0.002042000],\
                     [0.346200000,	0.396300000,	-2.582000000,	2.845000000,	-0.000280400],\
                     [2.883000000,	-5.873000000,	2.489000000,	1.510000000,	-0.002577000],\
                     [3.025000000,	-6.366000000,	3.137000000,	1.213000000,	-0.001367000],\
                     [3.229000000,	-6.844000000,	3.535000000,	1.088000000,	-0.002891000],\
                     [3.334000000,	-7.131000000,	3.829000000,	0.976600000,	-0.002952000],\
                     [3.146000000,	-6.855000000,	3.931000000,	0.786000000,	-0.002934000],\
                     [3.744000000,	-8.836000000,	6.018000000,	0.084070000,	0.000482500	]]
            
        RefCoef = [[16.320000,	-57.820000,	79.240000,	-50.080000,	13.340000],\
                     [40.480000,	-119.300000,	134.800000,	-70.970000,	16.110000],\
                     [57.490000,	-164.500000,	178.000000,	-88.750000,	18.840000],\
                     [5.714000,	-16.670000,	18.630000,	-9.756000,	3.074000],\
                     [-0.548800,	-6.498000,	21.200000,	-20.970000,	7.814000],\
                     [4.290000,	-12.670000,	14.660000,	-8.153000,	2.871000],\
                     [21.740000,	-64.440000,	74.890000,	-41.790000,	10.620000],\
                     [4.341000,	-12.800000,	14.780000,	-8.203000,	2.879000],\
                     [41.360000,	-117.800000,	127.600000,	-64.370000,	14.260000],\
                     [4.490000,	-12.660000,	13.970000,	-7.501000,	2.693000]]
        
        # Defining the 10 correlations 
        Ts_A = lambda Cos : TransCoef[0][0]*Cos**4+TransCoef[0][1]*Cos**3+TransCoef[0][2]*Cos**2+TransCoef[0][3]*Cos+TransCoef[0][4]
        Ts_B = lambda Cos : TransCoef[1][0]*Cos**4+TransCoef[1][1]*Cos**3+TransCoef[1][2]*Cos**2+TransCoef[1][3]*Cos+TransCoef[1][4]
        Ts_C = lambda Cos : TransCoef[2][0]*Cos**4+TransCoef[2][1]*Cos**3+TransCoef[2][2]*Cos**2+TransCoef[2][3]*Cos+TransCoef[2][4]
        Ts_D = lambda Cos : TransCoef[3][0]*Cos**4+TransCoef[3][1]*Cos**3+TransCoef[3][2]*Cos**2+TransCoef[3][3]*Cos+TransCoef[3][4]
        Ts_E = lambda Cos : TransCoef[4][0]*Cos**4+TransCoef[4][1]*Cos**3+TransCoef[4][2]*Cos**2+TransCoef[4][3]*Cos+TransCoef[4][4]
        Ts_F = lambda Cos : TransCoef[5][0]*Cos**4+TransCoef[5][1]*Cos**3+TransCoef[5][2]*Cos**2+TransCoef[5][3]*Cos+TransCoef[5][4]
        Ts_G = lambda Cos : TransCoef[6][0]*Cos**4+TransCoef[6][1]*Cos**3+TransCoef[6][2]*Cos**2+TransCoef[6][3]*Cos+TransCoef[6][4]
        Ts_H = lambda Cos : TransCoef[7][0]*Cos**4+TransCoef[7][1]*Cos**3+TransCoef[7][2]*Cos**2+TransCoef[7][3]*Cos+TransCoef[7][4]
        Ts_I = lambda Cos : TransCoef[8][0]*Cos**4+TransCoef[8][1]*Cos**3+TransCoef[8][2]*Cos**2+TransCoef[8][3]*Cos+TransCoef[8][4]
        Ts_J = lambda Cos : TransCoef[9][0]*Cos**4+TransCoef[9][1]*Cos**3+TransCoef[9][2]*Cos**2+TransCoef[9][3]*Cos+TransCoef[9][4]
        
        Rs_A = lambda Cos : RefCoef[0][0]*Cos**4+RefCoef[0][1]*Cos**3+RefCoef[0][2]*Cos**2+RefCoef[0][3]*Cos+RefCoef[0][4]
        Rs_B = lambda Cos : RefCoef[1][0]*Cos**4+RefCoef[1][1]*Cos**3+RefCoef[1][2]*Cos**2+RefCoef[1][3]*Cos+RefCoef[1][4]
        Rs_C = lambda Cos : RefCoef[2][0]*Cos**4+RefCoef[2][1]*Cos**3+RefCoef[2][2]*Cos**2+RefCoef[2][3]*Cos+RefCoef[2][4]
        Rs_D = lambda Cos : RefCoef[3][0]*Cos**4+RefCoef[3][1]*Cos**3+RefCoef[3][2]*Cos**2+RefCoef[3][3]*Cos+RefCoef[3][4]
        Rs_E = lambda Cos : RefCoef[4][0]*Cos**4+RefCoef[4][1]*Cos**3+RefCoef[4][2]*Cos**2+RefCoef[4][3]*Cos+RefCoef[4][4]
        Rs_F = lambda Cos : RefCoef[5][0]*Cos**4+RefCoef[5][1]*Cos**3+RefCoef[5][2]*Cos**2+RefCoef[5][3]*Cos+RefCoef[5][4]
        Rs_G = lambda Cos : RefCoef[6][0]*Cos**4+RefCoef[6][1]*Cos**3+RefCoef[6][2]*Cos**2+RefCoef[6][3]*Cos+RefCoef[6][4]
        Rs_H = lambda Cos : RefCoef[7][0]*Cos**4+RefCoef[7][1]*Cos**3+RefCoef[7][2]*Cos**2+RefCoef[7][3]*Cos+RefCoef[7][4]
        Rs_I = lambda Cos : RefCoef[8][0]*Cos**4+RefCoef[8][1]*Cos**3+RefCoef[8][2]*Cos**2+RefCoef[8][3]*Cos+RefCoef[8][4]
        Rs_J = lambda Cos : RefCoef[9][0]*Cos**4+RefCoef[9][1]*Cos**3+RefCoef[9][2]*Cos**2+RefCoef[9][3]*Cos+RefCoef[9][4]
    
        # Evaluating curves for each zone 
        
        alpha = np.linspace(10,90,9)
        Cosalpha = np.cos(np.deg2rad(alpha))
        
        if self.U <= 1.42 and self.SHGC > 0.45:                         #Zone 1
            Ts_alpha = Ts_E(Cosalpha)
            Rs_alpha = Rs_E(Cosalpha)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
        
        elif self.U <= 1.42 and self.SHGC <= 0.45  and self.SHGC > 0.35:      #Zone 2 linear interpolation
            Ts1 = Ts_J(Cosalpha)
            Rs1 = Rs_J(Cosalpha)
            Ts2 = Ts_E(Cosalpha)
            Rs2 = Rs_E(Cosalpha)            
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.SHGC - 0.35)/(0.45 - 0.35)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.SHGC - 0.35)/(0.45 - 0.35)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U <= 1.42 and self.SHGC <= 0.35   :              #Zone 3 
            Ts_alpha = Ts_J(Cosalpha)
            Rs_alpha = Rs_J(Cosalpha)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
           
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC > 0.55  :    #Zone 4 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = Ts_E(Cosalpha)
            Rs2 = Rs_E(Cosalpha)            
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 1.42)/(1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 1.42)/(1.7 - 1.42)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC <= 0.55  and self.SHGC > 0.5 :     #Zone 5 bilinear interpolation
            Ts11 = Ts_E(Cosalpha)
            Rs11 = Rs_E(Cosalpha)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs21 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts22 = Ts_E(Cosalpha)
            Rs22 = Rs_E(Cosalpha)
            Ts_alpha = (Ts11*(0.55- self.SHGC)*(1.7- self.U)+Ts21*( self.SHGC-0.5)*(1.7- self.U)+Ts12*(0.55- self.SHGC)*( self.U-1.42)+Ts22*( self.SHGC-0.5)*( self.U-1.42))/((0.55-0.5)*(1.7-1.42))
            Rs_alpha = (Rs11*(0.55- self.SHGC)*(1.7- self.U)+Rs21*( self.SHGC-0.5)*(1.7- self.U)+Rs12*(0.55- self.SHGC)*( self.U-1.42)+Rs22*( self.SHGC-0.5)*( self.U-1.42))/((0.55-0.5)*(1.7-1.42))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC <= 0.5  and self.SHGC > 0.45 :     #Zone 6 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)        
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 1.42)/(1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 1.42)/(1.7 - 1.42)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC <= 0.45  and self.SHGC > 0.35 :     #Zone 7 bilinear interpolation
            Ts11 = Ts_J(Cosalpha)
            Rs11 = Rs_J(Cosalpha)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs21 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts22 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs22 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts_alpha = (Ts11*(0.45- self.SHGC)*(1.7- self.U)+Ts21*( self.SHGC-0.35)*(1.7- self.U)+Ts12*(0.45- self.SHGC)*( self.U-1.42)+Ts22*( self.SHGC-0.35)*( self.U-1.42))/((0.45-0.35)*(1.7-1.42))
            Rs_alpha = (Rs11*(0.45- self.SHGC)*(1.7- self.U)+Rs21*( self.SHGC-0.35)*(1.7- self.U)+Rs12*(0.45- self.SHGC)*( self.U-1.42)+Rs22*( self.SHGC-0.35)*( self.U-1.42))/((0.45-0.35)*(1.7-1.42))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
    
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC <= 0.35  and self.SHGC > 0.3 :     #Zone 8 linear interpolation
            Ts1 = Ts_J(Cosalpha)
            Rs1 = Rs_J(Cosalpha)
            Ts2 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)         
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 1.42)/(1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 1.42)/(1.7 - 1.42)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC <= 0.3  and self.SHGC > 0.25 :     #Zone 9 bilinear interpolation
            Ts11 = Ts_J(Cosalpha)
            Rs11 = Rs_J(Cosalpha)
            Ts12 = Ts_J(Cosalpha)
            Rs12 = Rs_J(Cosalpha)
            Ts21 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)],axis=0)
            Rs21 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)],axis=0)
            Ts22 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs22 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts_alpha = (Ts11*(0.30- self.SHGC)*(1.7- self.U)+Ts21*( self.SHGC-0.25)*(1.7- self.U)+Ts12*(0.30- self.SHGC)*( self.U-1.42)+Ts22*( self.SHGC-0.25)*( self.U-1.42))/((0.30-0.25)*(1.7-1.42))
            Rs_alpha = (Rs11*(0.30- self.SHGC)*(1.7- self.U)+Rs21*( self.SHGC-0.25)*(1.7- self.U)+Rs12*(0.30- self.SHGC)*( self.U-1.42)+Rs22*( self.SHGC-0.25)*( self.U-1.42))/((0.30-0.25)*(1.7-1.42))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
        
        elif self.U > 1.42 and self.U <= 1.7 and self.SHGC <= 0.25  :                     #Zone 10 linear interpolation
            Ts1 = Ts_J(Cosalpha)
            Rs1 = Rs_J(Cosalpha)
            Ts2 = np.mean([Ts_F(Cosalpha),  Ts_H(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_F(Cosalpha),  Rs_H(Cosalpha)],axis=0)
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 1.42)/(1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 1.42)/(1.7 - 1.42)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U <= 3.41 and self.U > 1.7 and self.SHGC > 0.55  :                        #Zone 11
            Ts_alpha = Ts_E(Cosalpha)
            Rs_alpha = Rs_E(Cosalpha)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U <= 3.41 and self.U > 1.7 and self.SHGC <= 0.55 and self.SHGC > 0.5  :       #Zone 12
            Ts1 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs1 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts2 = Ts_E(Cosalpha)
            Rs2 = Rs_E(Cosalpha)                     
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.SHGC - 0.5)/(0.55 - 0.5)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.SHGC - 0.5)/(0.55 - 0.5)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U <= 3.41 and self.U > 1.7 and self.SHGC <= 0.5 and self.SHGC > 0.3   :     #Zone 13
            Ts_alpha = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs_alpha = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U <= 3.41 and self.U > 1.7 and self.SHGC <= 0.3 and self.SHGC > 0.25   :      #Zone 14             
            Ts1 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)],axis=0)
            Rs1 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)],axis=0)
            Ts2 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)          
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.SHGC - 0.25)/(0.3 - 0.25)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.SHGC - 0.25)/(0.3 - 0.25)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U <= 3.41 and self.U > 1.7 and self.SHGC <= 0.25    :                    #Zone 15
            Ts_alpha = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)],axis=0)
            Rs_alpha = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)],axis=0)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
                        
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC > 0.65     :                    #Zone 16 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = Ts_A(Cosalpha)
            Rs2 = Rs_A(Cosalpha)            
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 3.41)/(4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 3.41)/(4.54 - 3.41)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.65  and self.SHGC > 0.6  :    #Zone 17 bilinear interpolation
            Ts11 = Ts_E(Cosalpha)
            Rs11 = Rs_E(Cosalpha)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs21 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0) 
            Ts22 = Ts_A(Cosalpha)
            Rs22 = Rs_A(Cosalpha)
            Ts_alpha = (Ts11*(0.65- self.SHGC)*(4.54- self.U)+Ts21*( self.SHGC-0.6)*(4.54- self.U)+Ts12*(0.65- self.SHGC)*( self.U-3.41)+Ts22*( self.SHGC-0.6)*( self.U-3.41))/((0.65-0.6)*(4.54-3.41))
            Rs_alpha = (Rs11*(0.65- self.SHGC)*(4.54- self.U)+Rs21*( self.SHGC-0.6)*(4.54- self.U)+Rs12*(0.65- self.SHGC)*( self.U-3.41)+Rs22*( self.SHGC-0.6)*( self.U-3.41))/((0.65-0.6)*(4.54-3.41))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.6  and self.SHGC > 0.55:      #Zone 18 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 3.41)/(4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 3.41)/(4.54 - 3.41)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.55  and self.SHGC > 0.5 :     #Zone 19 bilinear interpolation
            Ts11 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs11 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs21 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts22 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs22 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts_alpha = (Ts11*(0.55- self.SHGC)*(4.54- self.U)+Ts21*( self.SHGC-0.5)*(4.54- self.U)+Ts12*(0.55- self.SHGC)*( self.U-3.41)+Ts22*( self.SHGC-0.5)*( self.U-3.41))/((0.55-0.5)*(4.54-3.41))
            Rs_alpha = (Rs11*(0.55- self.SHGC)*(4.54- self.U)+Rs21*( self.SHGC-0.5)*(4.54- self.U)+Rs12*(0.55- self.SHGC)*( self.U-3.41)+Rs22*( self.SHGC-0.5)*( self.U-3.41))/((0.55-0.5)*(4.54-3.41))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.5  and self.SHGC > 0.45  :    #Zone 20 linear interpolation
            Ts1 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs1 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts2 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 3.41)/(4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 3.41)/(4.54 - 3.41)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.45  and self.SHGC > 0.3  :    #Zone 21 bilinear interpolation
            Ts11 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs11 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts12 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs12 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts21 = Ts_D(Cosalpha)
            Rs21 = Rs_D(Cosalpha)
            Ts22 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs22 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts_alpha = (Ts11*(0.45- self.SHGC)*(4.54- self.U)+Ts21*( self.SHGC-0.3)*(4.54- self.U)+Ts12*(0.45- self.SHGC)*( self.U-3.41)+Ts22*( self.SHGC-0.3)*( self.U-3.41))/((0.45-0.3)*(4.54-3.41))
            Rs_alpha = (Rs11*(0.45- self.SHGC)*(4.54- self.U)+Rs21*( self.SHGC-0.3)*(4.54- self.U)+Rs12*(0.45- self.SHGC)*( self.U-3.41)+Rs22*( self.SHGC-0.3)*( self.U-3.41))/((0.45-0.3)*(4.54-3.41))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.3  and self.SHGC > 0.25 :     #Zone 22 bilinear interpolationm
            Ts11 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)],axis=0)
            Rs11 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)],axis=0)
            Ts12 = np.mean([Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)],axis=0)
            Rs12 = np.mean([Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)],axis=0)
            Ts21 = Ts_D(Cosalpha)
            Rs21 = Rs_D(Cosalpha)
            Ts22 = Ts_D(Cosalpha)
            Rs22 = Rs_D(Cosalpha)
            Ts_alpha = (Ts11*(0.3- self.SHGC)*(4.54- self.U)+Ts21*( self.SHGC-0.25)*(4.54- self.U)+Ts12*(0.3- self.SHGC)*( self.U-3.41)+Ts22*( self.SHGC-0.25)*( self.U-3.41))/((0.3-0.25)*(4.54-3.41))
            Rs_alpha = (Rs11*(0.3- self.SHGC)*(4.54- self.U)+Rs21*( self.SHGC-0.25)*(4.54- self.U)+Rs12*(0.3- self.SHGC)*( self.U-3.41)+Rs22*( self.SHGC-0.25)*( self.U-3.41))/((0.3-0.25)*(4.54-3.41))
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 3.41 and self.U <= 4.54 and self.SHGC <= 0.25 :                      #Zone 23 linear interpolation
            Ts1 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)],axis=0)
            Rs1 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)],axis=0)
            Ts2 = Ts_D(Cosalpha)
            Rs2 = Rs_D(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.U - 3.41)/(4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.U - 3.41)/(4.54 - 3.41)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 4.5  and self.SHGC > 0.65  :                                    #Zone 24            
            Ts_alpha = Ts_A(Cosalpha)
            Rs_alpha = Rs_A(Cosalpha)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 4.5 and self.SHGC > 0.6 and self.SHGC <= 0.65 :            #Zone 25 linear interpolation
            Ts1 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs1 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts2 = Ts_A(Cosalpha)
            Rs2 = Rs_A(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.SHGC - 0.6)/(0.65 - 0.6)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.SHGC - 0.6)/(0.65 - 0.6)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
        
        elif self.U > 4.5  and self.SHGC > 0.45 and self.SHGC <= 0.6:                    #Zone 26            
            Ts_alpha = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs_alpha = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
            
        elif self.U > 4.5  and self.SHGC > 0.3 and self.SHGC <= 0.45 :            #Zone 27 linear interpolation            
            Ts1 = Ts_D(Cosalpha)
            Rs1 = Rs_D(Cosalpha)
            Ts2 = np.mean([Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)],axis=0)
            Rs2 = np.mean([Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)],axis=0)
            Ts_alpha = Ts1 + (Ts2 - Ts1)*( self.SHGC - 0.3)/(0.45 - 0.3)
            Rs_alpha = Rs1 + (Rs2 - Rs1)*( self.SHGC - 0.3)/(0.45 - 0.3)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
    
        elif self.U > 4.5  and self.SHGC <= 0.3:                                 #Zone 28            
            Ts_alpha = Ts_D(Cosalpha)
            Rs_alpha = Rs_D(Cosalpha)
            Ts_abs_alpha = self.Ts*np.concatenate(([1],Ts_alpha))
            Rs_abs_alpha = self.Rs_f*np.concatenate(([1], Rs_alpha))
    
        else:
            print('Error in the Function. Check the Curve definition and areas')
            
        # Controlling the two vectors 
        
        Ts_abs_alpha[9] = 0
        Rs_abs_alpha_2 = np.zeros(10)
        Rs_abs_alpha_2[0] = Rs_abs_alpha[0]
        if max(Rs_abs_alpha) > 1:
            for i in range(1,10):
                Rs_abs_alpha_2[i] = Rs_abs_alpha[0] + (Rs_abs_alpha[i]-Rs_abs_alpha[0])*(1 - Rs_abs_alpha[0])/(max(Rs_abs_alpha)-Rs_abs_alpha[0])
            Rs_abs_alpha=Rs_abs_alpha_2
        else:
           Rs_abs_alpha[9]= 1  
        
        for i in range(10):
            if Rs_abs_alpha[i]+ Ts_abs_alpha[i] > 1:
                Rs_abs_alpha[i]= 1 - Ts_abs_alpha[i]   
                
        self.Rs_abs_alpha = Rs_abs_alpha
        self.Ts_abs_alpha = Ts_abs_alpha
        self.As_abs_alpha = 1 - self.Rs_abs_alpha - self.Ts_abs_alpha
        self.SHGC_abs_alpha = self.Ts_abs_alpha + self.N*self.As_abs_alpha
        self.alpha2 = np.linspace(0,90,10)
        self.SHGC_profile = interpolate.splrep(self.alpha2,self.SHGC_abs_alpha,s=0)

#%%--------------------------------------------------------------------------------------------------- 

'''
TEST METHOD
'''
import os

if __name__ =='__main__':
    
    Env=loadEvelopes(os.path.join('..','Input','buildings_envelope_V02_EP.xlsx'))
    
    '''
    x=np.linspace(0,90,10)
    import matplotlib.pyplot as plt
    for U in np.linspace(0.5,6.5,13):
        for SHGC in np.linspace(0.05,0.95,37):
            print('U='+str(U)+'\nSGHC='+str(SHGC))
            window=Window(pd.DataFrame({'name':'prova','U':U,'SHGC':SHGC,'tau_vis':0.9,'F_f':0.9,'F_sh':0.9},index=[0]))
            
            plt.plot(x,window.Rs_abs_alpha)
    '''