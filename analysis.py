import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.cluster import hierarchy
from tools.trees import dict_node
import math
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri, vectors
from rpy2.robjects.packages import STAP

def correlation(data, by_cols):
    '''
    Calculate a correlation matrix and returns a dict
    with the indication name the key and a matrix of dicts with
    Pearson r values P-values and count (the number of elements compared) as keys
    '''
    corr_array = []  
    if by_cols:
        x = np.hsplit(data.df_matrix, data.col_size)
    else:            
        x = np.vsplit(data.df_matrix, data.row_size)           
    for i in x:
        for j in x:
            mask = ~np.isnan(i) & ~np.isnan(j)
            if len(i[mask]) == 0:
                corr_array.append({
                    'pearson_r': None,
                    'p_value': None,
                    'count': None
                })
            else:
                slope, intercept, r_value, p_value, std_err = stats.linregress(i[mask],j[mask])
                corr_array.append({ 
                    'pearson_r': r_value,
                    'p_value': None if math.isnan(p_value) else p_value,
                    'count': len(i[mask])
                })

    #make the array into a matrix
    corr_array = np.array(corr_array)
    if by_cols:
        corr_mat = np.resize(corr_array,(data.col_size,data.col_size))
    else:
        corr_mat = np.resize(corr_array,(data.row_size,data.row_size))       
    # format output   
    response_output = {
        'name': data.name,
        'key': data.key,
        'matrix': {
          'rows': corr_mat.tolist(),
          'row_names': data.col_names if by_cols else data.row_names,
          'col_names': data.col_names if by_cols else data.row_names,
        }
    }
    return response_output

def z_score(data, by_cols):
    '''
    Calculate a z score matrix and returns a dict
    with the indication name the key and a matrix of dicts with
    z score values and count (the number of elements compared) as keys
    '''
    z_array = []  
    if by_cols:
        x = np.hsplit(data.df_matrix, data.col_size)
    else:            
        x = np.vsplit(data.df_matrix, data.row_size)           
    
    for arr in x:
        mask = ~np.isnan(arr)
        if len(arr[mask]) != 0:
            z_s =stats.zscore(arr[mask])
            arr = arr.tolist()
            for i, elem in enumerate(arr[0]):
                if math.isnan(elem):
                    z_s = np.insert(z_s,i,None,axis=None)
            z_array.append(z_s)    

    #make the array into a matrix
    z_array = np.array(z_array)
    if by_cols:
        z_mat = np.resize(z_array,(data.col_size,data.row_size))
    else:
        z_mat = np.resize(z_array,(data.row_size,data.col_size))
    z_mat = np.where(np.isnan(z_mat), None, z_mat)    
    # format output   
    response_output = {
        'name': data.name,
        'key': data.key,
        'matrix': {
          'rows': z_mat.tolist(),
          'row_names': data.row_names, 
          'col_names': data.col_names 
        }
    }
    return response_output       

def density(data,bandwidth):
    '''
    Uses gaussian kernel density estimation to return the  x_y values and
    count (number of items used) for a density plot 
    ''' 
    populations_array = np.vsplit(data.df_matrix,data.row_size) 
    density_array =[]
    for i,population in enumerate(populations_array):
        np.sort(population)
        mask = ~np.isnan(population)
        gaussian = gaussian_kde(population[mask],bw_method = bandwidth / population[mask].std(ddof=1))
        xs = np.linspace(0,1)
        density_array.append({
            'x_values': xs.tolist(),
            'density': gaussian(xs).tolist(),
            'count': len(population[mask])
        })
    #format output 
    response_output = {
        'name': data.name,
        'key': data.key,
        'matrix':{
          'row_names': data.row_names,
          'col_names': range(0,50),
          'rows': density_array
        }
    }
    return response_output    

def dendrogram(data,by_cols):
    dist_mat = data.to_distance_matrix(by_cols)
    clusters = hierarchy.average(dist_mat) 
    tree = hierarchy.to_tree(clusters, rd=False)
    leaf_labels =[]
    for i in leaves:
        if by_cols:
            leaf_labels.append(data.col_names[i])
        else:
            leaf_labels.append(data.row_names[i])
        
    
    response_output = {
        'name': data.name,
        'key': data.key,
        'tree': dict_node(tree, leaf_labels, 'root'),
        'labels': leaf_labels
    }
    return response_output

def DE(data,p_val,labels):
  
  pandas2ri.activate()
  # make labels into a dataframe -> tuple values ordered by dataframe column order
  labels_df = pd.DataFrame(labels)
  labels_df = labels_df.set_index('label')
  pd.to_numeric(labels_df['value'])
  df = data.df.T.join(labels_df)
  labels_tuple = tuple(list(df['value']))
  num = data.col_size*10
  rdf = pandas2ri.py2ri(data.df)
  with open('my_voom.R', 'r') as f:
    string = f.read()
  my_voom = STAP(string, "my_voom")
  iv= ro.IntVector(labels_tuple)
  group = ro.FactorVector(iv)
  fit = my_voom.run_voom(rdf,group,num)
  top = my_voom.topGenes (fit,2,pval=0.5)
  
  top_df = pandas2ri.ri2py(top)


  response_output = {
      'name': data.name,
      'key': data.key,
      'matrix':{
        'row_names': top_df.index.values.tolist(),
        'col_names': top_df.columns.values.tolist(),
        'rows': top_df.values.values.tolist()
      }
  }
  return response_output
  