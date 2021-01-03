from flask import Flask, render_template, request, redirect, url_for, jsonify, json
import flask

import numpy as np
import sys
import traceback
from request_handler import execute
from errors import ArchimedesError
from encoder import ArchimedesEncoder

app = Flask(__name__)

app.config.update(
    DEBUG=True,
    TRAP_HTTP_EXCEPTIONS=True,
    TRAP_BAD_REQUEST_ERRORS=True
)

@app.route('/')
def route_home():
    return "Archimedes is on"

@app.route('/', methods=['POST'])
def route_manifest():
    try:
        # post request formatted
        #post_request = request.get_json(force=True)
        data = execute(request)
        response = app.response_class(
            response=json.dumps(data, cls=ArchimedesEncoder),
            status=200,
            mimetype='application/json'
        )
        return response
        
    except Exception:
        print(traceback.format_exc())
        return flask.jsonify(msg = 500)
