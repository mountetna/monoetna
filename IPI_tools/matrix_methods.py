import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde
import pandas as pd
import math
from scipy.cluster import hierarchy
from scipy.spatial import distance
import json
from IPI_tools.tree_methods import add_node

class MatrixMethods:
    
    def __init__(self, name, key, row_names, col_names, rows,r_size,c_size):
        self.name = name
        self.key = key
        self.row_names = row_names
        self.col_names =col_names
        self.rows =rows
        self.r_size = r_size
        self.c_size = c_size
        return     
    
    def make_row_col_matrix(self):
        '''
        Extracts values matrices from request data
        '''''
        df = pd.DataFrame.from_items(zip(self.row_names,self.rows), orient='index', columns= self.col_names)
        series_matrix = df.as_matrix()
        series_matrix = np.array(series_matrix, dtype=np.float)
        return series_matrix
    

    def get_corr_matrix_data(self, by_cols):
        '''
        Calculate a correlation matrix and returns a dict
        with the indication name the key and a matrix of dicts with
        Pearson r values P-values and count (the number of elements compared) as keys
        '''
        series_matrix = self.make_row_col_matrix()
        corr_array = []  
        if by_cols:
            x = np.hsplit(series_matrix, self.c_size)
        else:            
            x = np.vsplit(series_matrix, self.r_size)           
        for i in x:
            for j in x:
                item_dict ={}
                mask = ~np.isnan(i) & ~np.isnan(j)
                if len(i[mask]) == 0:
                    item_dict['pearson_r'] = None
                    item_dict['p_value'] = None
                    item_dict['count'] = None
                    corr_array.append(item_dict)
                else:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(i[mask],j[mask])
                    item_dict['pearson_r'] = r_value            
                    if math.isnan(p_value): #returns none if P-value is NaN because jsonify could not handle NaNs
                        p_value = None
                    item_dict['p_value'] = p_value
                    item_dict['count'] = len(i[mask])
                    corr_array.append(item_dict)       
        #make the array into a matrix
        corr_array = np.array(corr_array)
        if by_cols:
            corr_mat = np.resize(corr_array,(self.c_size,self.c_size))
        else:
            corr_mat = np.resize(corr_array,(self.r_size,self.r_size))       
        # format output   
        response_output = {'name':self.name,'key':self.key,'matrix':{}}
        if by_cols:
            response_output['matrix']['row_names'] = self.col_names
            response_output['matrix']['col_names'] = self.col_names
        else:
            response_output['matrix']['row_names'] = self.row_names
            response_output['matrix']['col_names'] = self.row_names
        response_output['matrix']['rows'] =  corr_mat.tolist()
        return response_output         

    def get_density_plot_data(self,bandwidth):
        
        '''
        Uses gaussian kernel density estimation to return the  x_y values and
        count (number of items used) for a density plot 
        ''' 
        series = self.make_row_col_matrix()
        populations_array = np.vsplit(series,self.r_size) 
        density_array =[]
        for i,population in enumerate(populations_array):
            item_dict ={}
            np.sort(population)
            mask = ~np.isnan(population)
            density = gaussian_kde(population[mask],bw_method = bandwidth / population[mask].std(ddof=1))
            xs = np.linspace(0,1)
            item_dict['x_values'] = xs.tolist()
            item_dict['density'] = density(xs).tolist()
            item_dict['count'] = len(population[mask])
            density_array.append(item_dict)
        #format output 
        response_output = {'name': self.name,'key': self.key,'matrix':{}}
        response_output['matrix']['row_names'] = self.row_names
        response_output['matrix']['col_names'] = range(0,50)
        response_output['matrix']['rows'] =  density_array
        return response_output    
        
    def get_distance_matrix(self,by_cols):
        '''
        Calculates a symmetric distance matrix and returns the upper triangle
        '''
        series = self.make_row_col_matrix()
        dist_array = []
        if by_cols:
            x = np.hsplit(series,self.c_size)
        else:            
            x = np.vsplit(series,self.r_size)    
        for i in x:
            for j in x:
                mask = ~np.isnan(i) & ~np.isnan(j)
                if len(i[mask]) == 0:
                    dist_array.append(None)
                else:
                    dist_array.append(distance.euclidean(i[mask],j[mask]))

        dist_array = np.array(dist_array)
        if by_cols:
            dist_mat = np.resize(dist_array,(self.c_size,self.c_size))
        else:
            dist_mat = np.resize(dist_array,(self.r_size,self.r_size))
        # get the upper triangle of the symmetric distance matrix
        dist_mat = distance.squareform(dist_mat)         
        return dist_mat   
        
    def get_dendrogram_data(self,by_cols):
        
        dist_mat = self.get_distance_matrix(by_cols)
        clusters = hierarchy.average(dist_mat) 
        T = hierarchy.to_tree( clusters , rd=False )
        tree_dict = {'children':[],'name': 'root',}
        dendrogram = add_node(T,tree_dict,self.row_names)
        response_output = {'name': self.name,'key': self.key}
        response_output['dendrogram'] = dendrogram 
        return response_output
               
def main():
    
    print 'This is plot_method.py'
    

if __name__ == "__main__":
    main()

    