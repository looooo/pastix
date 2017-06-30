"""
 @file solver.py

 PaStiX python object oriented wrapper

 @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-06-29

"""
from .pastix import *
from .enum import *

class solver(object):

    def __init__(self, A=None, **kwargs):
        self.iparm, self.dparm = initParam()
        self.verbose = kwargs.setdefault("verbose", 1)
        if not self.verbose: # 0 or False
            self.iparm[iparm.verbose] = verbose.Not
        elif self.verbose==2:
            self.iparm[iparm.verbose] = verbose.Yes
        else: # 1 or True or anything else, default value
            self.iparm[iparm.verbose] = verbose.No
        self.pastix_data = init( self.iparm, self.dparm )
        if A is not None:
            self.setup(A)

    def setup(self, A):
        """
        Setup the solver object for the classic case w/o Schur complement
        """
        self.A = A
        self.spmA = spm(A)
        if self.verbose:
            self.spmA.printInfo()
        task_analyze(self.pastix_data, self.spmA)
        task_numfact(self.pastix_data, self.spmA)

    def solve(self, b, refine=True, x0=None, check=False):
        """
        Solve and refine the full problem Ax = b with refinement when no Schur Complement is involved
        """
        x = b.copy()
        task_solve(self.pastix_data, x)
        if refine:
            task_refine(self.pastix_data, b, x)
        if check:
            self.spmA.checkAxb(x0, b, x)
        return x

    def schur(self, A, schur_list):
        """
        Setup the solver object to work in Schur complement mode
        """

        self.factotype = factotype.LLT
        self.iparm[iparm.factorization] = self.factotype
        self.A    = A
        self.spmA = spm(A)

        if self.verbose:
            self.spmA.printInfo()

        schur_list = np.asarray(schur_list, pastix_int)
        self.schur_list = schur_list + self.spmA.findBase()
        setSchurUnknownList(self.pastix_data, schur_list)

        task_analyze(self.pastix_data, self.spmA)
        task_numfact(self.pastix_data, self.spmA)

        nschur = len(schur_list)
        self.S = np.zeros( (nschur, nschur), order='F', dtype=A.dtype )
        getSchur( self.pastix_data, self.S )

    def schur_forward(self, b):
        """
        Solves the forward step of the Schur problem:
               A x = b with A = L S L^h

        This step solves  L f = b with f = S L^h x
        """
        x = b.copy()
        # 1- Apply P to b
        subtask_applyorder( self.pastix_data, direction.Forward, x )

        # 2- Forward solve on the non Schur complement part of the system
        if self.factotype == factotype.LU:
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

    def schur_backward(self, y, b, refine=True, x0=None, check=False):
        """
        Solves the backward step of the Schur problem:
               A x = b with A = L S L^h

        This step solves L^h x = y with y the solution of L S y = b
        """
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

        if check:
            self.spmA.checkAxb(x0, b, x)
        return x

    def finalize(self):
        finalize(self.pastix_data, self.iparm, self.dparm)
