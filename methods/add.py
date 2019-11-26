from helper import helper

@helper("Sum of all the input arguments.")
def func(*args):
    return sum(args)