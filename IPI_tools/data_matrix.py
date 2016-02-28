from IPI_tools.matrix_methods import MatrixMethods

class DataMatrix:
    def __init__(self,series):

        self.name = series['name']
        self.key = series['key']
        self.row_names = series['matrix']['row_names']
        self.col_names = series['matrix']['col_names']
        self.rows = series['matrix']['rows']
        return
    

    def matrix(self):
        '''
        Extracts rows, row names ,column names request data to
        make a PlottingMethods object
        '''
        
        data = MatrixMethods(name = self.name, key = self.key,
                                   row_names = self.row_names,
                                   col_names = self.col_names,
                                   rows = self.rows, r_size = len(self.row_names),
                                   c_size = len(self.col_names))
        return data

def main():
    print 'This is parse_json_input.py'

if __name__ == "__main__":
    main()    