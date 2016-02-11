import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import pandas as pd
import math


def get_data(input_request):
    '''
    Extracts ans array of values matrices from request data
    '''''
    series_array = []
    for s in input_request['input']['series']:
        matrix = s['matrix']
        df = pd.DataFrame.from_items(zip(matrix['row_names'],matrix['rows']), orient='index', columns= matrix['col_names'])
        series_matrix = df.as_matrix()
        series_matrix = np.array(series_matrix, dtype=np.float)
        series_array.append(series_matrix)
    return series_array

def get_corr_matrix(input_request):
    
    '''
    Calculate a correlation matrix from input request
    returns an array of matrices of dicts with R_squared values P-values and
    the number of elements compared as keys
    '''

    series_array = get_data(input_request)
    by_cols = input_request['params']['columns']
    corr_mat_array =[]
    for sample_array in series_array:
        columns = np.size(sample_array,1)
        rows = np.size(sample_array,0)
        corr_array = []

        if by_cols:
            x = np.hsplit(sample_array,columns)
        else:            
            x = np.vsplit(sample_array,rows)    
        
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
                    #returns none if P-value is NaN because jsonify could not handle NaNs
                    if math.isnan(p_value):
                        p_value = None
                    item_dict['p_value'] = p_value
                    item_dict['count'] = len(i[mask])
                    corr_array.append(item_dict)
        
        corr_array = np.array(corr_array)
        if by_cols:
            corr_mat = np.resize(corr_array,(columns,columns))
        else:
            corr_mat = np.resize(corr_array,(rows,rows))    
    corr_mat_array.append(corr_mat)
    
    # formatting output to look like the request   
    new_input ={} 
    for i, s in enumerate(input_request['input']['series']):
        new_input['series'] = []
        new_input['series'].append({'name':s['name'],'key':s['key'],'matrix':{}})
        if by_cols:
            new_input['series'][i]['matrix']['col_names'] = s['matrix']['col_names']
            new_input['series'][i]['matrix']['row_names'] = s['matrix']['col_names']
        else:
            new_input['series'][i]['matrix']['col_names'] = s['matrix']['row_names']
            new_input['series'][i]['matrix']['row_names'] = s['matrix']['row_names']
        new_input['series'][i]['matrix']['rows'] =  corr_mat_array[i].tolist()
    return new_input

    
def main():
    print 'This is plot_method.py'
    

if __name__ == "__main__":
    main()

    