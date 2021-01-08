from pandas import Series, DataFrame
import functions_vector
import functions_matrix

def add_to_library(module):
    namespace = vars(module)
    public = (name for name in namespace if name[:1] != "_")
    for name in getattr(module, "__all__", public):
        setattr(Library, name, namespace[name])

class LibraryFunction:
    def __init__(self, func):
        self.func = func

    def call(self, *args):
        return self.func(*args)

#   [ :diff_exp, :pca, :center, :scale, :transpose, :sd, :correlation, :wilcox, :beeswarm, :normal, :density ]

class Library:
    def help():
        return 'Help'

add_to_library(functions_vector)
add_to_library(functions_matrix)

