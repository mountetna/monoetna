from flask import Flask,render_template, request, redirect, url_for, jsonify, json
import flask
import urllib
from IPI_tools import compensation
from IPI_tools.compenstaion import get_PCA
from IPI_tools.plot_methods import get_corr_matrix
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
        if input_request['params']['method'] == 'pca':
            pca_mat = np.array(input['matrix'])
            pca_fracs, pca_Wt = get_PCA(pca_mat)
            return flask.jsonify(var_vec=pca_fracs.tolist(),vectors=pca_Wt.tolist())
        if input_request['params']['method'] == 'correlation':
            corr_mat = get_corr_matrix(input_request)
            return flask.jsonify(corr_mat=corr_mat)                 
        return flask.jsonify(error='Method not found:'+request.form['json_text']), 422

    except Exception:
        type, value, tb = sys.exc_info()
        return flask.jsonify(type=str(type),
                             error=str(value),
                             traceback=traceback.extract_tb(tb))

