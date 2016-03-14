import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.cluster import hierarchy
from tools.trees import dict_node
import math

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
    if by_cols:
        labels = data.col_names
    else:
        labels = data.row_names
    response_output = {
        'name': data.name,
        'key': data.key,
        'tree': dict_node(tree, labels, 'root')
    }
    return response_output
