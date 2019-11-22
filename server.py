from flask import Flask, render_template, request, redirect, url_for, jsonify, json
import flask

import numpy as np
import sys
import traceback
import importlib
import methods

app = Flask(__name__)
app.config.update(
    DEBUG=True,
    TRAP_HTTP_EXCEPTIONS=True,
    TRAP_BAD_REQUEST_ERRORS=True
)

@app.route('/')
def route_test():
    return "test page"


"""
input JSON string example
{
  "func": "function_name",
  "args": [
    1, 2, 3, 4
  ]
}
"""

@app.route('/json/', methods=['POST'])
def route_json():
    
    try:
        # post request formatted
        post_request = request.get_json(force=True)
        func_name = post_request["func"]
        
        args = []
        if "args" in post_request:
            if isinstance(post_request["args"], type(None)):
                pass
            else:    
                args = post_request["args"]
            
            
        method = importlib.import_module("%s.%s"%("methods",func_name))
        return flask.jsonify(method.func(*args))
            
    except Exception:
        type, value, tb = sys.exc_info()
        return flask.jsonify(type=str(type),
                         error=str(value),
                         traceback=traceback.extract_tb(tb))

