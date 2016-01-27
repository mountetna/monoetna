import FlowCytometryTools
from FlowCytometryTools import FCMeasurement
import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from matplotlib import pyplot as plt

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

def get_r_squared_corr_matrix(sample_array, cols):
    
    r_squared =[]
    columns = np.size(sample_array,1)
    rows = np.size(sample_array,0)
    if cols:       
        x = np.hsplit(sample_array,columns)
        for i in x:
            for j in x:
                slope, intercept, r_value, p_value, std_err = stats.linregress(i.flat,j.flat)
                r_squared.append(r_value ** 2)
        r_squared = np.array(r_squared)
        r_squared = np.resize(r_squared,(columns,columns))
        return r_squared
    else:            
        x = np.vsplit(sample_array,rows)
        for i in x:
            for j in x:
                slope, intercept, r_value, p_value, std_err = stats.linregress(i.flat,j.flat)
                r_squared.append(r_value ** 2)
        r_squared = np.array(r_squared)
        r_squared = np.resize(r_squared,(rows,rows))
        return r_squared
    
    
    
def get_eigenvalues_and_eigenvectors(sample_array):
    
    cov_mat = np.cov(sample_array, rowvar=0)
    return np.linalg.eig(cov_mat)
    
def main():
    
    datafile= fcs_file
    sample = FCMeasurement(ID='', datafile=datafile)
    get_PCA_FCS(sample)
    

    

if __name__ == "__main__":
    main()

    