import numpy as np
import pandas as pd
from scipy.spatial import distance
import math
from sklearn.preprocessing import Imputer

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
    
    def filter_zero_rows_cols(self, how='All'):
      '''
      will remove either rows with all zeroes or columns with all zeroes
      or both
      '''
      self.df = self.df.replace(0,np.nan)
      if how == 'rows':
        self.df = self.df.dropna(how='all',axis=0)
      elif how == 'columns': 
        self.df = self.df.dropna(how='all',axis=1)
      else:
        self.df = self.df.dropna(how='all',axis=0)
        self.df = self.df.dropna(how='all',axis=1)
      return self.df.replace(np.nan,0)
        
    def to_distance_matrix(self,by_cols):
      '''
      Calculates a euclidean distance matrix and returns the upper triangle
      '''
      # compress the matrix to remove rows and cols that are all null
      comp_df = self.df.dropna(axis=(0), how='all')
      comp_df_matrix = np.array(comp_df.as_matrix(),dtype=np.float)
      # replace nulls with median of the array
      imp = Imputer(missing_values='NaN',strategy="median",axis=1)
      imp.fit(comp_df_matrix)
      imputed_matrix = imp.transform(comp_df_matrix)
      if by_cols:
        return distance.pdist(imputed_matrix.transpose(),metric='correlation')
      else:    
        return distance.pdist(imputed_matrix,metric='correlation')  
