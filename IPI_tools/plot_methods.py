import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde
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
    count (the number of elements compared) as keys
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
    
    # format output   
    response_output ={} 
    for i, s in enumerate(input_request['input']['series']):
        response_output['series'] = []
        response_output['series'].append({'name':s['name'],'key':s['key'],'matrix':{}})
        if by_cols:
            response_output['series'][i]['matrix']['row_names'] = s['matrix']['col_names']
            response_output['series'][i]['matrix']['col_names'] = s['matrix']['col_names']
        else:
            response_output['series'][i]['matrix']['row_names'] = s['matrix']['row_names']
            response_output['series'][i]['matrix']['col_names'] = s['matrix']['row_names']
        response_output['series'][i]['matrix']['rows'] =  corr_mat_array[i].tolist()
    return response_output  
    
    

def get_density_plot_data(input_request):
    
    '''
    Uses gaussian kernel density estimation to return the  x_y values and
    count (number of items used) for a density plot 
    ''' 
    series_array = get_data(input_request)
    bandwidth = input_request['params']['bandwidth']
    density_plots_array = []
    for i, series in enumerate(series_array):
        row_names = input_request['input']['series'][i]['matrix']['row_names']
        rows = np.size(series,0)
        populations_array = np.vsplit(series,rows) 
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

    density_plots_array.append(density_array)    
    
    response_output ={} 
    for i, s in enumerate(input_request['input']['series']):
        response_output['series'] = []
        response_output['series'].append({'name':s['name'],'key':s['key'],'matrix':{}})
        response_output['series'][i]['matrix']['row_names'] = s['matrix']['row_names']
        response_output['series'][i]['matrix']['col_names'] = range(0,50)
        response_output['series'][i]['matrix']['rows'] =  density_plots_array[i]
    return response_output    

          
    
def main():
    print 'This is plot_method.py'
    

if __name__ == "__main__":
    main()

    