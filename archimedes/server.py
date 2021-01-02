from flask import Flask, render_template, request, redirect, url_for, jsonify, json
import flask

import numpy as np
import sys
import traceback
from request_handler import execute
from errors import ArchimedesError

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
        return flask.jsonify(execute(request))
        
    except Exception:
        print(traceback.format_exc())
        return flask.jsonify(msg = 500)
