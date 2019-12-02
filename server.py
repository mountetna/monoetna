from flask import Flask, render_template, request, redirect, url_for, jsonify, json
import flask

import numpy as np
import sys
import traceback
from request_handler import execute
from errors import FunctionError, TetherError

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
        return flask.jsonify(status = 200, 
                             output = execute(post_request))
        
    except FunctionError as e:
        return flask.jsonify(msg = e.msg,
                             status = e.status,
                             info = e.info)
    
    except TetherError as e:
        return flask.jsonify(msg = e.msg,
                             status = e.status,
                             info = e.info,
                             request_info = e.request_info)
    
    except Exception:
        print(traceback.format_exc())
        return flask.jsonify(msg = 500)
