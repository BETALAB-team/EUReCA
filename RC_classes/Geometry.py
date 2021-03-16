'''IMPORTING MODULES'''

import os
import numpy as np
import pyclipper as pc
from RC_classes.Envelope import loadEnvelopes
from RC_classes.auxiliary_functions import wrn

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions

def unit_normal(a, b, c):
    
    '''
    This function starting from three points defines the normal vector of the plane
    
    Parameters
    ----------
    a : list of floats
        list with three floats with the coordinates of the first point
    b : list of floats
        list with three floats with the coordinates of the second point
    c : list of floats
        list with three floats with the coordinates of the third point
 
    Returns
    -------
    tuple: the three components of the normal vector (anticlockwise)
    
    '''
    
    # Check input data type
    
    for ipt in (a,b,c):
        if not isinstance(ipt, list) and not isinstance(ipt, tuple):
            raise TypeError(f'ERROR unit_normal function, an input is not a list: input {ipt}')
        if len(ipt) != 3:
            raise TypeError(f'ERROR unit_normal function, a vertex is not a list of 3 components: input {ipt}')
        try:
            ipt = list(ipt)
            ipt[0] = float(ipt[0])
            ipt[1] = float(ipt[1])
            ipt[2] = float(ipt[2])
        except ValueError:
            raise ValueError(f'ERROR unit_normal function, a coordinate is not a float: input {ipt}')
    
    # unit normal vector of plane defined by points a, b, and c
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#%%

def poly_area(poly):
   
    '''
    From a list of points calculates the area of the polygon
    
    Parameters
    ----------
    poly : list of list of floats (polygon)
        list with three floats with the coordinates of the first point
 
    Returns
    -------
    float: the polygon area
    
    '''
    
    # Check input data type
    
    if not isinstance(poly, list):
        raise TypeError(f'ERROR poly_area function, the input is not a list: input {poly}')
    for vtx in poly:
        if not isinstance(vtx, list) and not isinstance(vtx, tuple):
            print(type(vtx))
            raise TypeError(f'ERROR poly_area function, an input is not a list: input {vtx}')
        if len(vtx) != 3:
            raise TypeError(f'ERROR poly_area function, a vertex is not a list of 3 components: input {vtx}')
        try:
            vtx = list(vtx)
            vtx[0] = float(vtx[0])
            vtx[1] = float(vtx[1])
            vtx[2] = float(vtx[2])
        except ValueError:
            raise ValueError(f'ERROR unit_normal function, a coordinate is not a float: input {vtx}')

    # area of polygon poly
    
    if len(poly) < 3:                                                          # Not a plane - no area
        wrn('WARNING number of vertices lower than 3, the area will be zero')
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return float(abs(result/2))

#%%

def check_complanarity(vertListTot,precision=1):
    '''
    checks the complanarity of a list of points
    
    Parameters
    ----------
    vertListTot : list of list of floats (polygon)
        list with n lists of three floats (n vertices)
    precision: float
        defines the precision of the control
 
    Returns
    -------
    boolean: are they in the same plane? True or False 
    
    '''
    
    # Check input data type
    
    if not isinstance(vertListTot, list):
        raise TypeError(f'ERROR check_complanarity function, the input is not a list: input {vertListTot}')
    for vtx in vertListTot:
        if not isinstance(vtx, list):
            raise TypeError(f'ERROR check_complanarity function, an input is not a list: input {vtx}')
        if len(vtx) != 3:
            raise TypeError(f'ERROR check_complanarity function, a vertex is not a list of 3 components: input {vtx}')
        try:
            vtx[0] = float(vtx[0])
            vtx[1] = float(vtx[1])
            vtx[2] = float(vtx[2])
        except ValueError:
            raise ValueError(f'ERROR check_complanarity function, a coordinate is not a float: input {vtx}')    
    
    try:
        precision = float(precision)
    except ValueError:
            raise ValueError(f'ERROR check_complanarity function, precision is not a float: precision {precision}')    
    
    # Look if they are coplanar
    
    flag = True
    for i in range(len(vertListTot)-3):
        vertList = vertListTot[i:(i+4)]
        a1 = vertList[1][0] - vertList[0][0]
        b1 = vertList[1][1] - vertList[0][1]
        c1 = vertList[1][2] - vertList[0][2]
        a2 = vertList[2][0] - vertList[0][0]
        b2 = vertList[2][1] - vertList[0][1]
        c2 = vertList[2][2] - vertList[0][2]
        a = b1 * c2 - b2 * c1
        b = a2 * c1 - a1 * c2
        c = a1 * b2 - b1 * a2
        d = (- a * vertList[0][0] - b * vertList[0][1] - c * vertList[0][2])

        # equation of plane is: a*x + b*y + c*z = 0
        # checking if the 4th point satisfies
        # the above equation
        if not (np.abs(a * vertList[3][0] + b * vertList[3][1] + c * vertList[3][2] + d) < precision):
            flag = False

    return flag

#%%

def centroid(vertList):
    
    '''
    From a list of points calculates the centroid
    
    Parameters
    ----------
    vertList : list of list of floats (polygon)
        list with n lists of three floats (n vertices)
 
    Returns
    -------
    np.array: coordinates of the centroid (3 components)
    
    '''
    
    # Check input data type
    
    if not isinstance(vertList, list):
        raise TypeError(f'ERROR centroid function, the input is not a list: input {vertList}')
    for vtx in vertList:
        if not isinstance(vtx, list):
            raise TypeError(f'ERROR centroid function, an input is not a list: input {vtx}')
        if len(vtx) != 3:
            raise TypeError(f'ERROR centroid function, a vertex is not a list of 3 components: input {vtx}')
        try:
            vtx[0] = float(vtx[0])
            vtx[1] = float(vtx[1])
            vtx[2] = float(vtx[2])
        except ValueError:
            raise ValueError(f'ERROR centroid function, a coordinate is not a float: input {vtx}')    
      
      
    # Centroid calculation
    c = np.array([0.,0,0], dtype = float)
    for i in vertList:
        c += np.array(i)
    return c/len(vertList)

#%%

def normalAlternative(vertList):
    
    '''
    Alternative
    This function starting from three points defines the normal vector of the plane
    
    Parameters
    ----------
    vertList : list of list of floats (polygon)
        list with n lists of three floats (n vertices)
 
    Returns
    -------
    np.array: coordinates of the centroid (3 components)
    
    '''
    
    # Check input data type
    
    if not isinstance(vertList, list):
        raise TypeError(f'ERROR normalAlternative function, the input is not a list: input {normalAlternative}')
    # if len(vertList) > 3:
        # wrn(f"normalAlternative function, there vertlist should be 3 components long: vertList {vertList}")
        
    for vtx in vertList:
        if not isinstance(vtx, list):
            raise TypeError(f'ERROR normalAlternative function, an input is not a list: input {vtx}')
        if len(vtx) != 3:
            raise TypeError(f'ERROR normalAlternative function, a vertex is not a list of 3 components: input {vtx}')
        try:
            vtx[0] = float(vtx[0])
            vtx[1] = float(vtx[1])
            vtx[2] = float(vtx[2])
        except ValueError:
            raise ValueError(f'ERROR normalAlternative function, a coordinate is not a float: input {vtx}')    
      
      
    c = centroid(vertList)
    crossProd=np.array([0.,0,0])
    for i in range(len(vertList)):
        a = np.array(vertList[i-1]) - c
        b = np.array(vertList[i]) - c
        crossProd += np.cross(a,b)
    return crossProd / np.linalg.norm(crossProd)

#%%

def project(x,proj_axis):
    
    '''
    Internal Function
      
    # Project onto either the xy, yz, or xz plane. (We choose the one that avoids degenerate configurations, which is the purpose of proj_axis.)
    # In this example, we would be projecting onto the xz plane.
    '''
    return tuple(c for i, c in enumerate(x) if i != proj_axis)

#%%

def project_inv(x,proj_axis,a,v):
    
    '''
    Internal Function
      
    # Returns the vector w in the walls' plane such that project(w) equals x.
    '''
    w = list(x)
    w[proj_axis:proj_axis] = [0.0]
    c = a
    for i in range(3):
        c -= w[i] * v[i]
    c /= v[proj_axis]
    w[proj_axis] = c
    return tuple(w)

#%%--------------------------------------------------------------------------------------------------- 
#%% Surface class

class Surface:
    
    '''
    Class surface checks the complanarity and calculates the area.
    Then calculates the azimuth and tilt of the surface and set a surface
    type depending on the tilt angle
    
    Methods:
    __init__:
        name
        sky dome azimuth subdivision
        sky dome height subdivision
        window wall ratio
        list of the vertices
        surface type: ExtWall, Roof, GroundFloor, Ceiling, IntWall    
        
    maxHeight: no input
    
    minHeight: no input
    
    checkSurfaceCoincidence takes a second surface to check if they are co-planar and they match eachother:
        secondary surface
        
    reduceArea reduce the area of the surface:
        area to reduce
        
    Methods:
        __init__
        maxHeight
        minHeight
        checkSurfaceCoincidence
        reduceArea
    
    '''
    
    def __init__(self,name,azSubdiv,hSubdiv,wwr,rh_gross,vertList=[[0,0,0],[0,0,0],[0,0,0]],surfType='ExtWall'):
        
        '''
        input:
            list of vertices: [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],..]
            typology: one of these strings: 'ExtWall','Roof','GroundFloor'

        complanarity:
            https://www.geeksforgeeks.org/program-to-check-whether-4-points-in-a-3-d-plane-are-coplanar/

        the area is calculated from:
            https://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates
        
        Parameters
            ----------
            name : string
                name of the surface
            azSubdiv: int
                number of azimuth subdivision
            hSubdiv: int
                number of tilt subdivision
            wwr : list of floats
                four floats with the window to wall ratio (N E S W)
            rh_gross : float
                multiplicator of the surface area
            vertList : list of list of floats (polygon)
                list with n lists of three floats (n vertices)
            surfType : string
                string that defines the surface type.
                'ExtWall' or 'GroundFloor' or 'Roof' 
                
        Returns
        -------
        None.        
        
        '''

        # Check input data type
        
        if not isinstance(name, str):
            raise TypeError(f'ERROR Surface class geometry, name is not a string: name {name}') 
        if not isinstance(azSubdiv, int):
            raise TypeError(f'ERROR Surface class geometry, azSubdiv is not a int: azSubdiv {azSubdiv}') 
        if not isinstance(hSubdiv, int):
            raise TypeError(f'ERROR Surface class geometry, azSubdiv is not a int: azSubdiv {hSubdiv}') 
        if not isinstance(wwr, list) or not isinstance(wwr[0], float) or not isinstance(wwr[1], float) or not isinstance(wwr[2], float) or not isinstance(wwr[3], float):
            raise TypeError(f'ERROR Surface class geometry, wwr is not a list of floats: wwr {wwr}') 
        if not isinstance(rh_gross, float):
            raise TypeError(f'ERROR Surface class geometry, rh_gross is not a float: rh_gross {rh_gross}') 
        if not surfType in ['ExtWall' , 'GroundFloor' , 'Roof' ]:
            raise TypeError(f'ERROR Surface class geometry, surfType is not a correct string: surfType {surfType}')
            
        if not isinstance(vertList, list):
            raise TypeError(f'ERROR Surface class geometry , the input is not a list: input {vertList}')
        for vtx in vertList:
            if not isinstance(vtx, list):
                raise TypeError(f'ERROR Surface class geometry, an input is not a list: input {vtx}')
            if len(vtx) != 3:
                raise TypeError(f'ERROR Surface class geometry, a vertex is not a list of 3 components: input {vtx}')
            try:
                vtx[0] = float(vtx[0])
                vtx[1] = float(vtx[1])
                vtx[2] = float(vtx[2])
            except ValueError:
                raise ValueError(f'ERROR Surface class geometry, a coordinate is not a float: input {vtx}')          

        # Check input data quality
        
        if azSubdiv > 10 or hSubdiv > 5:
            wrn(f"WARNING Surface class, init, solar calculation could be long..... azSubdiv {azSubdiv}, hSubdiv {hSubdiv}")
        for wwr_ in wwr:
            if wwr_ < 0. or wwr_ > 0.9:
                wrn(f"WARNING Surface class, init, are you sure about the window to wall ratio ? wwr {wwr}")
        if rh_gross < 0.5 or rh_gross > 1.5:
            wrn(f"WARNING Surface class, init, are you sure about the external walls multiplication coeff?? rh_gross {rh_gross}")
           
        self.name = name

        # Check coplanarity

        if not check_complanarity(vertList):
            print('surface points not planar')
            self.type = '--'
            self.area = 0.
            self.normal = 0.
            self.height = 0.
            self.azimuth = 0.
            return
        
        self.OnOff_shading = 'Off'
        self.type = surfType
        self.vertList = vertList
        
        # Area calculation
        
        self.area = poly_area(self.vertList)
        if self.area == 0.:
            self.area = 0.0000001
        
        '''
        Considering only three points in calculating the normal vector could create
        reverse orientations if the three points are in a non-convex angle of the surface

        for this reason theres an alternative way to calculate the normal,
        implemented in function: normalAlternative

        reference: https://stackoverflow.com/questions/32274127/how-to-efficiently-determine-the-normal-to-a-polygon-in-3d-space
        '''
        
        #self.normal = unit_normal(self.vertList[0], self.vertList[1], self.vertList[2])
        self.normal = normalAlternative(self.vertList)

        # set the azimuth and zenith

        if self.normal[2] == 1:
            self.height = 0
            self.azimuth = 0
        elif self.normal[2] == -1:
            self.height= 180
            self.azimuth = 0
        else:
            self.height = 90 - np.degrees(np.arctan((self.normal[2]/(np.sqrt(self.normal[0]**2+self.normal[1]**2)))))
            if self.normal[1] == 0:
                if self.normal[0] > 0:
                    self.azimuth = -90
                elif self.normal[0] < 0:
                    self.azimuth = 90
            else:
                if self.normal[1]<0:
                    self.azimuth = np.degrees(np.arctan(self.normal[0]/self.normal[1]))
                else:
                    if self.normal[0] < 0:
                        self.azimuth = 180 + np.degrees(np.arctan(self.normal[0]/self.normal[1]))
                    else:
                        self.azimuth = -180 + np.degrees(np.arctan(self.normal[0]/self.normal[1]))
        
        # Set surface inclination
        
        if self.height < 40:
            self.type = 'Roof'
        if self.height > 150:
            self.type = 'GroundFloor'
        if self.type == 'ExtWall':
            self.area = self.area*rh_gross
        
        # Azimuth and tilt approximation 
        
        delta_a = 360/(2*azSubdiv)
        delta_h = 90/(2*hSubdiv)
        x = np.arange(-delta_h,90+2*delta_h,2*delta_h)
        
        for n in range(len(x)-1):
            if self.height >= x[n] and self.height < x[n+1]:
                self.height_round = int((x[n]+x[n+1])/2)
                self.F_r = (1+np.cos(np.radians(self.height_round)))/2                                                          
            elif self.height >= x[-1] and self.height < 150:
                self.height_round = 90
                self.F_r = (1+np.cos(np.radians(self.height_round)))/2
            else:
                self.height_round = 0                                          # Only to avoid errors           
                
        y = np.arange(-180-delta_a,180+2*delta_a,2*delta_a)
        for n in range(len(y)-1):
            if self.azimuth >= y[n] and self.azimuth < y[n+1]:
                self.azimuth_round = int((y[n]+y[n+1])/2)
                if self.azimuth_round == 180:
                    self.azimuth_round = -180
                    
        if self.height_round == 0:
            self.azimuth_round = 0
        
        # Set the window area
        
        if self.type == 'ExtWall':
            self.centroid_coord = centroid(vertList)
            if 135 < self.azimuth_round <= 180 or -180 <= self.azimuth_round < -135:
                self.wwr = wwr[0]
            elif -135 <= self.azimuth_round <= -45:
                self.wwr = wwr[1]
            elif -45 < self.azimuth_round < 45:
                self.wwr = wwr[2]
            elif 45 <= self.azimuth_round <= 135:
                self.wwr = wwr[3]
        else:
            self.wwr = 0
            
        self.opaqueArea = (1-self.wwr)*self.area
        self.glazedArea = (self.wwr)*self.area
        if self.glazedArea == 0:
            self.glazedArea = 0.0000001                                        #Avoid zero division
        if self.opaqueArea == 0:
            self.opaqueArea = 0.0000001                                        #Avoid zero division
            
                        
    def maxHeight(self):
    
        '''
        Find the higher vertex

        
        Parameters
            ----------
            None
                
        Returns
        -------
        float.    [m]     
        
        '''
    
        hmax = 0
        for vert in self.vertList:
            hmax = max(hmax,vert[2])
        return hmax


    def minHeight(self):
    
        '''
        Find the lower vertex

        
        Parameters
            ----------
            None
                
        Returns
        -------
        float.        [m] 
        
        '''
    
        hmin = 10000
        for vert in self.vertList:
            hmin = min(hmin,vert[2])
        return hmin


    def checkSurfaceCoincidence(self,otherSurface):
        """
        Check if two surface are coincident    
        
        Parameters
            ----------
            otherSurface : EUReCA.RC_classes.geometry.Surface
                another surface object
                
        Returns
        -------
        boolean. Are the surfaces coincident? True/Flase         
        """
        
        # Check Input data type
        
        if not isinstance(otherSurface,Surface):
            raise ValueError(f"ERROR Surface class, surface {self.name}, checkSurfaceCoincidence. otherSurface is not a Surface object: otherSurface {otherSurface}")
        
        # Check the coincidence of two surface looking firstly at the coplanarity 
        # of the points and then the direction of the normals vectors
        
        flagPoints = False
        plane = self.vertList
        
        # Coplanarity test
        for i in otherSurface.vertList:
            if check_complanarity(plane + [i],precision=5):
                flagPoints = True
        
        # Normal vector test
        flagNormal = False
        if np.linalg.norm(self.normal+otherSurface.normal)<0.2:
            flagNormal= True
        return (flagNormal and flagPoints)

    def calculateIntersectionArea(self,otherSurface):
        '''
        Claculates the area between two adjacent surfaces
        
        reference: https://stackoverflow.com/questions/39003450/transform-3d-polygon-to-2d-perform-clipping-and-transform-back-to-3d
        
        Parameters
            ----------
            otherSurface : EUReCA.RC_classes.geometry.Surface
                another surface object
                
        Returns
        -------
        float. The area [m2]
        '''
        
        # Check Input data type
        
        if not isinstance(otherSurface,Surface):
            raise ValueError(f"ERROR Surface class, surface {self.name}, calculateIntersectionArea. otherSurface is not a Surface object: otherSurface {otherSurface}")
        
        # Check the coincidence of two surface looking firstly at the coplanarity 
        # of the points and then the direction of the normals vectors
        a = self.normal[0]*self.vertList[0][0]+self.normal[1]*self.vertList[0][1]+self.normal[2]*self.vertList[0][2]
        proj_axis = max(range(3), key=lambda i: abs(self.normal[i]))
        projA = [project(x,proj_axis) for x in self.vertList]
        projB = [project(x,proj_axis) for x in otherSurface.vertList]
        scaledA = pc.scale_to_clipper(projA)
        scaledB = pc.scale_to_clipper(projB)
        clipper = pc.Pyclipper()
        clipper.AddPath(scaledA, poly_type=pc.PT_SUBJECT, closed=True)
        clipper.AddPath(scaledB, poly_type=pc.PT_CLIP, closed=True)
        intersections = clipper.Execute(pc.CT_INTERSECTION, pc.PFT_NONZERO, pc.PFT_NONZERO)
        intersections = [pc.scale_from_clipper(i) for i in intersections]
        if len(intersections)==0:
            return 0
        intersection = [project_inv(x,proj_axis,a,self.normal) for x in intersections[0]]
        area = poly_area(intersection)
        return area if area>0 else 0.


    def reduceArea(self,AreaToReduce):
        '''
        Reduces the area of the surface
        
        Parameters
            ----------
            AreaToReduce : float
                the area to subtract [m2]
                
        Returns
        -------
        None.   
        '''
        
        # Check Input data type
        
        if not isinstance(AreaToReduce, float) or AreaToReduce < 0.:
            try:
                AreaToReduce = float(AreaToReduce)
            except ValueError:    
                raise ValueError(f"ERROR Surface class, surface {self.name}, reduceArea. The area is not a positive float: AreaToReduce {AreaToReduce}")
        
        # Area reduction
        
        if self.area - AreaToReduce > 0.0000001:
            self.area = self.area - AreaToReduce
            self.opaqueArea = (1-self.wwr)*self.area
            self.glazedArea = (self.wwr)*self.area
        else:
            self.area = 0.0000001
            self.opaqueArea = 0.0000001
            self.glazedArea = 0.0000001


    def printInfo(self):
    
        '''
        Just prints some attributes of the surface
        
        Parameters
            ----------
            None
                
        Returns
        -------
        None.   
        '''
    
        print('Name: '+self.name+\
              '\nArea: '+str(self.area)+\
              '\nType: '+str(self.type)+\
              '\nAzimuth: '+str(self.azimuth)+\
              '\nHeight: '+str(self.height)+\
              '\nVertices: '+str(self.vertList)+\
              '\n')


#%%--------------------------------------------------------------------------------------------------- 
#%% SurfaceInternalMass class

class SurfaceInternalMass():
    
    '''
    Class to define a surface for thermal capacity using area and surface type
    with a specific geometry
    
    Methods:
        init
    '''

    def __init__(self,name,area=0.0000001,surfType='IntWall'):
        '''
        input:
            area: area of the internal surface
            surfType: 'IntWall' or 'IntCeiling'

        attrubutes:
            area
            surfType

        Parameters
            ----------
            name : string
                name of the surface
            area: float
                number of azimuth subdivision [m2]
            surfType : string
                string that defines the surface type.
                'IntWall' or  'IntCeiling'  or 'IntFloor'
                
        Returns
        -------
        None.        
        
        '''

        # Check input data type
        
        if not isinstance(name, str):
            raise TypeError(f'ERROR SurfaceInternalMass class geometry, name is not a string: name {name}') 
        if not isinstance(area, float) or area < 0.:
            raise TypeError(f'ERROR SurfaceInternalMass class geometry, area is not a positive float: area {area}') 
        if not surfType in ['IntWall' , 'IntCeiling' ,'IntFloor']:
            raise TypeError(f'ERROR SurfaceInternalMass class geometry, surfType is not a correct string: surfType {surfType}')

        # Sets some attributes
        
        self.name = name
        self.area = area
        self.type = surfType
        if self.area < 0.0000001:
            self.area = 0.0000001
        self.opaqueArea = self.area

#%%--------------------------------------------------------------------------------------------------- 
#%% SurfaceInternalAdjacent class

class SurfaceInternalAdjacent(SurfaceInternalMass):
    
    '''
    Inherited from SurfaceInternalMass
    adds the adjacentZone attribute
    Currently not used in the code
    
    '''
    def __init__(self,name,area,surfType='IntCeiling',adjacentZone = None):
        super().__init__(area,surfType)
        self.adjacentZone = adjacentZone

#%%--------------------------------------------------------------------------------------------------- 
#%% some test methods

'''
TEST METHOD
'''
'''
if __name__ == '__main__':
    surf=Surface([[0,0,0],[1,1,0],[2,0,1],[4,0,2],[2,-2,2]])
    print('area '+str(surf.area))
    print('normal '+str(surf.normal))
    print('height '+str(surf.height))
    print('azimuth '+str(surf.azimuth))
    print(surf.glazedArea)
'''
'''
if __name__ == '__main__':
    surf=Surface('surf1',[[726604.119, 5032389.067, 65.808],
                  [726611.079, 5032388.257, 65.808],
                  [726610.939, 5032387.267, 65.808],
                  [726613.499, 5032386.907, 65.808],
                  [726613.649, 5032387.977, 65.808],
                  [726626.989, 5032386.107, 65.808],
                  [726626.869, 5032385.227, 65.808],
                  [726629.498, 5032384.857, 65.808],
                  [726629.628, 5032385.787, 65.808],
                  [726636.528, 5032384.818, 65.808],
                  [726637.768, 5032393.677, 65.808],
                  [726605.129, 5032397.697, 65.808]])

    print('normal '+str(surf.normal))
    print('alternative normal '+str(surf.normal))
    surf1=Surface('surf2',[[726604.119, 5032389.067, 55.018],
                 [726605.129, 5032397.697, 55.018],
                 [726637.768, 5032393.677, 55.018],
                 [726636.528, 5032384.817, 55.018],
                 [726629.628, 5032385.787, 55.018],
                 [726629.498, 5032384.857, 55.018],
                 [726626.868, 5032385.227, 55.018],
                 [726626.988, 5032386.107, 55.018],
                 [726613.649, 5032387.977, 55.018],
                 [726613.499, 5032386.907, 55.018],
                 [726610.939, 5032387.267, 55.018],
                 [726611.079, 5032388.257, 55.018]])
    print('normal '+str(surf1.normal))
    print('alternative normal '+str(surf1.normal))

    print(surf1.checkSurfaceCoincidence(surf))

'''
'''
    surf=SurfaceInternalMass(45,'IntWall')
    print(surf.area)
    print(surf.surfType)

    surf=SurfaceInternalAdjacent(45,'IntWall','bsjaf')
    print(surf.area)
    print(surf.surfType)
    print(surf.adjacentZone)
'''
'''
if __name__ == '__main__':
    surf=Surface('surf1',[[2,1,0],[4,3,0],[3,2,2],[2,1,2]])

    #print('normal '+str(surf.normal))
    #print('alternative normal '+str(surf.normal))


    surf1=Surface('surf2',[[4,2,0],[4,2,1],[6,4,1],[6,4,0]])
    #print('normal '+str(surf1.normal))
    #print('alternative normal '+str(surf1.normal))

    #print(surf1.checkSurfaceCoincidence(surf))
    
    env_path = os.path.join('..','Input','buildings_envelope_V02_EP.xlsx')
    envelopes = loadEvelopes(env_path)
'''
if __name__ == '__main__':
    surf=Surface('Surf1',[[655338.4496301126, 5033359.473662503, 6.0], [655338.4496301126, 5033359.473662503, 6.0], [655347.8363386589, 5033352.433631093, 6.0], [655352.9990283594, 5033358.691436791, 6.0]])
    
