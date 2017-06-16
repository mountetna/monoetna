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
    series = input_request['input']['series'][0] # assuming only 1 series
    data = DataMatrix(series)
    pythia_output ={}
    for key, parameters in params.iteritems():
      method = parameters['method']

      if method == 'correlation':
        by_cols = parameters['columns']
        corr_mat = {
          'series': [
            analysis.correlation(data,by_cols)
          ]
        }
        pythia_output[key] = corr_mat
      
      elif method == 'density':
        bandwidth =parameters['bandwidth']
        density_plots = {
          'series': [
            analysis.density(data,bandwidth)
          ]
        }
        pythia_output[key]=density_plots
      
      elif method == 'dendrogram':
        by_cols = parameters['columns']
        dendrogram = analysis.dendrogram(data,by_cols)
        pythia_output[key]=dendrogram
    
      elif method == 'z_score':
        by_cols = parameters['columns']
        z_mat = {
          'series': [
            analysis.z_score(data,by_cols)
          ]
        }
        pythia_output[key] = z_mat
          
      elif method == 'DE':
        p_val = parameters['p_value']
        labels = parameters['labels']
        DE_table = {
          'series': [
            analysis.DE(data,p_val,labels)
          ]
        }
        pythia_output[key] = DE_table    
      else:       
        return flask.jsonify(error='Method not found:'+request.form['json_text']), 422

    return flask.jsonify(pythia_output)
  except Exception:
    type, value, tb = sys.exc_info()
    return flask.jsonify(type=str(type),
                         error=str(value),
                         traceback=traceback.extract_tb(tb))

