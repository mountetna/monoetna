from flask import Flask,render_template, request, redirect, url_for, jsonify, json
import flask
import urllib
from tools.data_matrix import DataMatrix

import analysis
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
        params = input_request['params']
        method = params['method']
        for i, series in enumerate(input_request['input']['series']):

            data = DataMatrix(series)

            if method == 'correlation':
                by_cols = params['columns']
                corr_mat = {
                  'series': [
                    analysis.correlation(data,by_cols)
                  ]
                }
                return flask.jsonify(corr_mat=corr_mat)
            
            elif method == 'density':
                bandwidth =params['bandwidth']
                density_plots = {
                  'series': [
                    analysis.density(data,bandwidth)
                  ]
                }
                return flask.jsonify(density_plots=density_plots)
            
            elif method == 'dendrogram':
                by_cols = params['columns']
                dendrogram = analysis.dendrogram(data,by_cols)
                return flask.jsonify(dendrogram=dendrogram)
                
        
        return flask.jsonify(error='Method not found:'+request.form['json_text']), 422

    except Exception:
        type, value, tb = sys.exc_info()
        return flask.jsonify(type=str(type),
                             error=str(value),
                             traceback=traceback.extract_tb(tb))

