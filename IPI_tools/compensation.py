import FlowCytometryTools
from FlowCytometryTools import FCMeasurement
import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from matplotlib import pyplot as plt
from scipy import stats
import pandas as pd
import math

def get_inverse_spillover_matrix(sample):
    
    spillover = sample.meta['SPILL']
    spillover = spillover.split(',')
    spillover = np.array(spillover)
    para_num = int(spillover[0])
    spillover = np.delete(spillover,[0]) # remove 1st element which is the number of parameters in the matrix 
    spillover = np.reshape(spillover,(para_num+1,para_num))
    len_para_headers = len(spillover[0])
    spillover = np.delete(spillover,0, axis=0)# remove the headers
    return (len_para_headers, np.linalg.inv(spillover.astype(float)))

def get_compensated_array(sample):
    
    len_para_headers,spillover_inv = get_inverse_spillover_matrix(sample)
    sample_array = sample.data.values
    channel_names = sample.channel_names
    len_sample_headers = len(sample.channel_names)
    # split vertically to obtain columns with parameters that need to be compensated
    split_sample_array = np.split(sample_array,[len_sample_headers-len_para_headers], axis=1)
    #the dot product of the sample sub array and the spillover matrix = compensated parameter values
    split_comp_array = np.dot(split_sample_array[1],spillover_inv)
    compensated_sample_array = np.concatenate((split_sample_array[0],split_comp_array),axis =1)
    hlog_comp_array = FlowCytometryTools.core.transforms.hlog(compensated_sample_array)
    return hlog_comp_array
    

def get_PCA(sample_array):
    
    PCA = mlabPCA(sample_array)
    return (PCA.fracs, PCA.Wt)
    
def get_PCA_FCS(temp_file_path):
    
    datafile= temp_file_path
    sample = FCMeasurement(ID='', datafile=datafile)
    hlog_comp_array = get_compensated_array(sample)
    pca_fracs,pca_Wt = get_PCA(hlog_comp_array)
    return (pca_fracs,pca_Wt)

def get_correlation_data(input_request):
    '''
    Extracts values matrix from request data
    '''''
    input_params = input_request['input']['series'][0] 
    matrix = input_params['matrix']
    df = pd.DataFrame.from_items(zip(matrix['row_names'],matrix['rows']), orient='index', columns= matrix['col_names'])
    sample_array = df.as_matrix()
    sample_array = np.array(sample_array, dtype=np.float)
    return sample_array
    

def get_corr_matrix(input_request):
    
    '''
    Calculate a correlation matrix from input request
    returns a matrix of dicts with R_squared values P-values and
    the number of elements compared as keys
    '''

    sample_array = get_correlation_data(input_request)
    by_cols = input_request['params']['columns']
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
                item_dict['r_squared'] = None
                item_dict['p_value'] = None
                item_dict['num_items_compared'] = None
                corr_array.append(item_dict)
            else:
                slope, intercept, r_value, p_value, std_err = stats.linregress(i[mask],j[mask])
                item_dict['r_squared'] = r_value**2
                #returns none if P-value is NaN because jsonify could not handle NaNs
                if math.isnan(p_value):
                    p_value = None
                item_dict['p_value'] = p_value
                item_dict['num_items_compared'] = len(i[mask])
                corr_array.append(item_dict)
    
    corr_array = np.array(corr_array)
    if by_cols:
        corr_mat = np.resize(corr_array,(columns,columns))
    else:
        corr_mat = np.resize(corr_array,(rows,rows))    

    # formatting output to look like the request   
    new_input ={} 
    for i, s in enumerate(input_request['input']['series']):
        new_input['series'] = []
        new_input['series'].append({'name':s['name'],'key':s['key'],'matrix':{}})
        new_input['series'][i]['matrix']['col_names'] = s['matrix']['col_names']
        new_input['series'][i]['matrix']['row_names'] = s['matrix']['row_names']
        new_input['series'][i]['matrix']['rows'] =  corr_mat.tolist()
    return new_input

  
    
def get_eigenvalues_and_eigenvectors(sample_array):
    
    cov_mat = np.cov(sample_array, rowvar=0)
    return np.linalg.eig(cov_mat)
    
def main():
    
    datafile= fcs_file
    sample = FCMeasurement(ID='', datafile=datafile)
    get_PCA_FCS(sample)
    

if __name__ == "__main__":
    main()

    