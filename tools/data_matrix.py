import numpy as np
import pandas as pd
from scipy.spatial import distance
import math

class DataMatrix:
    def __init__(self,series):
        self.name = series['name']
        self.key = series['key']
        self.row_names = series['matrix']['row_names']
        self.col_names = series['matrix']['col_names']
        self.rows = series['matrix']['rows']
        self.row_size = len(self.row_names)
        self.col_size = len(self.col_names)

        df = pd.DataFrame.from_items(zip(self.row_names,self.rows), 
            orient='index', columns= self.col_names)
        self.df = np.array(df.as_matrix(), dtype=np.float)
        
    def to_distance_matrix(self,by_cols):
        '''
        Calculates a symmetric distance matrix and returns the upper triangle
        '''
        dist_array = []
        if by_cols:
            x = np.hsplit(self.df,self.col_size)
        else:            
            x = np.vsplit(self.df,self.row_size)    
        for i in x:
            for j in x:
                mask = ~np.isnan(i) & ~np.isnan(j)
                if len(i[mask]) == 0:
                    dist_array.append(None)
                else:
                    dist_array.append(distance.euclidean(i[mask],j[mask]))

        dist_array = np.array(dist_array)

        if by_cols:
            dist_mat = np.resize(dist_array,(self.col_size,self.col_size))
        else:
            dist_mat = np.resize(dist_array,(self.row_size,self.row_size))

        # get the upper triangle of the symmetric distance matrix
        dist_mat = distance.squareform(dist_mat)
        return dist_mat   
