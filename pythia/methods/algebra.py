import numpy as np
from helper import helper

@helper("Addition of numbers")
def add(a, b):
    return a + b

@helper("Subtraction of numbers")
def sub(a, b):
    return a - b
    
@helper("Multiplication")
def mult(a, b):
    return a * b
    
@helper("Division")
def div(a, b):
    return a / b

@helper("Exponentiation")
def exp(a, b):
    return a ** b

@helper("Logarithm")
def log(a, b = np.e):
    return np.log(a)/np.log(b)