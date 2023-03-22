# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 12:18:34 2020

@author: pratenr15640
"""
import os
import sys
import json
import pandas as pd

class Cityjson_obj():
    def __init__(self,file):
        with open(file, 'r') as f:
            self.city = json.load(f)
        
    def addAttribute(self,dictionary,object_type,attribute):
        for obj_key in self.city['CityObjects'].keys():
            obj = self.city['CityObjects'][obj_key]
            if obj['type']== object_type:
                try:
                    obj['attributes'][attribute] = dictionary[obj_key]
                except KeyError:
                    sys.exit('Cityjson keys and dictionary keys are not the same')
            else:
                continue
            
    def saveJson(self,path):
        with open(path,'w') as file:
            json.dump(self.city,file)
            
    def reduceBdNumber(self,bd_list):
        new_obj = {}
        for k in bd_list:
            new_obj[k] = self.city['CityObjects'][k]
        print(new_obj)
        self.city['CityObjects'] = new_obj
        
        
        

if __name__=='__main__':
    p = os.path.join('DatiCityJSON','PaduaClipped_reprojected.json')
    data = pd.read_excel('dati_statistici edifici_v2.xlsx', sheet_name = 'ListaPiovego', usecols = 'B:F', skiprows = 2, header = 1)
    labels_uso = pd.Series(data['Tipologia Immobile'].values,index=data.Edifici).to_dict()
    labels_età = pd.Series(data['Age'].values,index=data.Edifici).to_dict()

    city =  Cityjson_obj(p)
        
    city.addAttribute(labels_uso,'Building','Use')
    city.addAttribute(labels_età,'Building','Age')
    city.saveJson(os.path.join('DatiCityJSON','Belzoni3.json'))
    
    