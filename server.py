from flask import Flask, render_template, request, redirect, url_for, jsonify, json
import flask

import numpy as np
import sys
import traceback
from request_handler import execute

app = Flask(__name__)
app.config.update(
    DEBUG=True,
    TRAP_HTTP_EXCEPTIONS=True,
    TRAP_BAD_REQUEST_ERRORS=True
)

@app.route('/')
def route_home():
    return "Home"

@app.route('/json/', methods=['POST'])
def route_json():
    
    try:
        # post request formatted
        post_request = request.get_json(force=True)
        return flask.jsonify(execute(post_request))
            
    except Exception:
        dtype, value, tb = sys.exc_info()
        return flask.jsonify(type=str(dtype),
                         error=str(value),
                         traceback=traceback.format_exc())

