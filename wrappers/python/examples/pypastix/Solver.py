from .pastix import *

class Solver(object):

    def __init__(self, A=None, **kwargs):
        self.parameters = kwargs
        self.iparm, self.dparm = initParam()
        self.pastix_data = init( self.iparm, self.dparm )
        if A is not None:
            self.setup(A)

    def setup(self, A):
        self.A = A
        self.spmA = spm(A)
        if self.parameters.setdefault("verbose", False):
            self.spmA.printInfo()
        analyze(self.pastix_data, self.spmA)
        numfact(self.pastix_data, self.spmA)

    def solve(self, b, refine_=True, x_check=None):
        x = b.copy()
        solve(self.pastix_data, self.spmA, x)
        if refine_:
            refine(self.pastix_data, b, x)
        if x_check is not None:
            self.spmA.checkAxb(x_check, b, x)
        return x

    def schur(self, A, schur_list):
        self.factotype = factotype.LLT
        self.iparm[iparm.factorization] = self.factotype
        self.A = A
        self.spmA = spm(A)
        if self.parameters.setdefault("verbose", False):
            self.spmA.printInfo()
        schur_list = np.asarray(schur_list, pastix_int)
        self.schur_list = schur_list
        setSchurUnknownList(self.pastix_data, schur_list)
        analyze(self.pastix_data, self.spmA)
        numfact(self.pastix_data, self.spmA)
        nschur = len(schur_list)
        self.S = np.zeros( (nschur, nschur), order='F', dtype=A.dtype )
        getSchur( self.pastix_data, self.S )

    def b2f(self, b):
        x = b.copy()
        # 1- Apply P to b
        subtask_applyorder( self.pastix_data, direction.Forward, x )

        # 2- Forward solve on the non Schur complement part of the system
        if 0 and self.factotype == factotype.LU:
            subtask_trsm( self.pastix_data, side.Left,
                                 uplo.Lower, trans.NoTrans,
                                 diag.Unit, x )
        else:
            subtask_trsm( self.pastix_data, side.Left,
                                 uplo.Lower, trans.NoTrans,
                                 diag.NonUnit, x )

        self.x = x
        nschur = len(self.schur_list)
        return x[-nschur:]

    def y2x(self, y, b, refine_=True, x_check=None):
        nschur = len(self.schur_list)
        x = self.x
        x[-nschur:] = y
        # 4- Backward solve on the non Schur complement part of the system
        if self.factotype == factotype.LU:
            subtask_trsm( self.pastix_data, side.Left,
                                 uplo.Upper, trans.NoTrans,
                                 diag.NonUnit, x )
        else:
            subtask_trsm( self.pastix_data, side.Left,
                                 uplo.Lower, trans.ConjTrans,
                                 diag.NonUnit, x )

        #  5- Apply P^t to x
        subtask_applyorder( self.pastix_data, direction.Backward, x )
        if refine_:
            refine(self.pastix_data, b, x)
        if x_check is not None:
            self.spmA.checkAxb(x_check, b, x)
        return x

    def finalize(self):
        finalize(self.pastix_data, self.iparm, self.dparm)
