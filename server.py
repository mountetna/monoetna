from flask import Flask,render_template, request, redirect, url_for, jsonify, json
import flask
import urllib
from IPI_tools import compensation
from IPI_tools.data_matrix import DataMatrix
from IPI_tools.matrix_methods import MatrixMethods
from IPI_tools.compensation import get_PCA
import numpy as np
import os
import sys
import traceback

app = Flask(__name__)
app.config.update(
    DEBUG=True,
    TRAP_HTTP_EXCEPTIONS=True,
    TRAP_BAD_REQUEST_ERRORS=True
)
@app.route('/test/')
def route_test():
    return flask.render_template('template.html')

@app.route('/json/', methods=['POST'])
def route_json():
    
    try:
        input_request = request.get_json(force=True)
        method = input_request['params']['method']
        for i, input in enumerate(input_request['input']['series']):
            data = DataMatrix(input).matrix()
            if method == 'pca':
                pca_mat = np.array(input['matrix'])
                pca_fracs, pca_Wt = get_PCA(pca_mat)
                return flask.jsonify(var_vec=pca_fracs.tolist(),vectors=pca_Wt.tolist())
            
            elif method == 'correlation':
                by_cols = input_request['params']['columns']
                corr_mat = {'series': []}
                corr_mat['series'].append(data.get_corr_matrix_data(by_cols))
                return flask.jsonify(corr_mat=corr_mat)
            
            elif method == 'density':
                bandwidth =input_request['params']['bandwidth']
                density_plots = {'series': []}
                density_plots['series'].append(data.get_density_plot_data(bandwidth))
                return flask.jsonify(density_plots=density_plots)
            
            elif method == 'dendrogram':
                by_cols = input_request['params']['columns']
                dendrogram = {'series': []}
                dendrogram['series'].append(data.make_dendrogram_data(by_cols))
                return flask.jsonify(dendrogram=dendrogram)
                
        
        return flask.jsonify(error='Method not found:'+request.form['json_text']), 422

    except Exception:
        type, value, tb = sys.exc_info()
        return flask.jsonify(type=str(type),
                             error=str(value),
                             traceback=traceback.extract_tb(tb))

