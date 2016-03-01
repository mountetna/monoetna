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

        self.df = pd.DataFrame.from_items(zip(self.row_names,self.rows), 
            orient='index', columns= self.col_names)
        self.df_matrix = np.array(self.df.as_matrix(), dtype=np.float)
        
    def to_distance_matrix(self,by_cols):
        '''
        Calculates a symmetric distance matrix and returns the upper triangle
        '''
        dist_array = []
        # compress the matrix to remove rows and cols that are all null
        comp_df = self.df.dropna(axis=(0,1), how='all')
        comp_df_matrix = np.array(comp_df.as_matrix(),dtype=np.float)
        comp_row_size = comp_df.shape[0]
        comp_col_size = comp_df.shape[1]
        # replace nulls with median of the array 
        imputed_matrix =[]
        for row in comp_df_matrix:
            median = np.median(row[~np.isnan(row)])
            row[np.isnan(row)] = median
            imputed_matrix.append(row)
        imputed_matrix = np.array(imputed_matrix)

        if by_cols:
            x = np.hsplit(imputed_matrix,comp_col_size)
        else:            
            x = np.vsplit(imputed_matrix,comp_row_size)    
        for i in x:
            for j in x:   
                dist_array.append(distance.euclidean(i,j))
        dist_array = np.array(dist_array)
       
        if by_cols:
            dist_mat = np.resize(dist_array,(comp_col_size,comp_col_size))
        else:
            dist_mat = np.resize(dist_array,(comp_row_size,comp_row_size))
        # get the upper triangle of the symmetric distance matrix
        dist_mat = distance.squareform(dist_mat)
        return dist_mat 
