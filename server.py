from flask import Flask,render_template, request, redirect, url_for, jsonify, json
import flask
import urllib
from IPI_tools import compensation
from IPI_tools.compensation import get_PCA_FCS, get_PCA, get_r_squared_corr_matrix
import os
import numpy as np

app = Flask(__name__)

@app.route('/test/')
def route_test():
    return flask.render_template('template.html')

@app.route('/json/', methods=['POST'])
def route_json():
    params = json.loads(request.form['json_text'])
    
    #if params['method'] == 'pca-FCS':
    #    FCSfile = urllib.URLopener()
    #    FCSfile.retreive("http://localhost:5000/json/blah", "temp.fcs")
    #    path = os.getcwd()
    #    pca_fracs, pca_Wt = get_PCA_FCS(path+'/temp.fcs')
    #    os.remove(path+'/temp.fcs')
    #    return flask.jsonify(var_vec=pca_fracs, vectors= pca_Wt)
    
    if params['method'] == 'pca':
        pca_mat = np.array(params['matrix'])
        pca_fracs, pca_Wt = get_PCA(pca_mat)
        return flask.jsonify(var_vec=pca_fracs.tolist(),vectors=pca_Wt.tolist())
    if params['method'] == 'correlation':
        corr_mat = np.array(params['matrix'])
        cols = params['columns']
        r_squared_mat = get_r_squared_corr_matrix(corr_mat, cols)
              
        
    return flask.jsonify(error='Method not found:'+request.form['json_text']), 422
    
if __name__ == '__main__':
    app.debug = True
    app.run(host='0.0.0.0', port=5000)
