
class ARBMatrixCerealWrap:
    """
    A wrapper class to enable serialization of complex arb matrix objects. The original arb matrix
    can be constructed via the `arb_matrix` method.
    """
    def __init__(self, arb_mat):
        self.nrows = arb_mat.nrows()
        self.ncols = arb_mat.ncols()
        self.arb_entries = [ [x.mid(), x.diameter()] for x in arb_mat.list()]
        self.base_ring = arb_mat.base_ring()

    def ball_field_elem(self, x):
        return self.field(x[0]).add_error(x[1])
    
    def arb_matrix(self):
        return matrix(self.base_ring, self.nrows , self.ncols,
                      map(self.ball_field_elem, self.arb_entries)  )
        
    def list(self):
        return self.arb_entries

