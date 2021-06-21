#
# This file allows us to generate:
#      - Documentation files
#      - IPARM/DPARM and their related enums declaration file
#        ( $PASTIX_HOME/include/pastix/api.h )
#      - parse_iparm, parse_dparm and parse_enums implementation and declaration.
#        pastix_ENUM_getstr implementation and declaration.
#        ( $PASTIX_HOME/common/parse_options.[h/c] )W
#
# If you want to modify one of these files, please modify this one.
#
###
# SYNTAX documentation. [] are used for optional keys.
#
# enum_name = { Name of the enum. Will be declared as pastix_[enum_name]_t.
#   "doc" : {
#     "brief": short description of the enum.
#     ["details"]: Elaborate the description. Used for the documentation.(1)
#   }
#   "values" : {
#     Value count [int] : {
#       "name": Enum string (Must begin with Pastix...)
#       ["value"]: To define the value of the enum [int]
#       ["brief"]: Short description of the impact of this value.
#     }
#   }
# }
###

enums = []

task = {
	"name" : "task",
	"doc" : {
		"brief" : "Main steps for the pastix() interface.",
		"details" : r'''
Those enums are used of the IPARM_START_TASK and IPARM_END_TASK parameters
that configure the pastix() call.
''',
	},
	"values" : {
		0 : {
    		"name": "PastixTaskInit",
    		"value": "0",
    		"brief": "Startup the library"
    	},
		1 : {
    		"name": "PastixTaskOrdering",
    		"value": "1",
    		"brief": "Ordering"
    	},
		2 : {
    		"name": "PastixTaskSymbfact",
    		"value": "2",
    		"brief": "Symbolic factorization"
    	},
		3 : {
    		"name": "PastixTaskAnalyze",
    		"value": "3",
    		"brief": "Tasks mapping and scheduling"
    	},
		4 : {
    		"name": "PastixTaskNumfact",
    		"value": "4",
    		"brief": "Numerical factorization"
    	},
		5 : {
    		"name": "PastixTaskSolve",
    		"value": "5",
    		"brief": "Numerical solve"
    	},
		6 : {
    		"name": "PastixTaskRefine",
    		"value": "6",
    		"brief": "Numerical refinement"
    	},
		7 : {
    		"name": "PastixTaskClean",
    		"value": "7",
    		"brief": "Clean"
    	}
	}
}

enums.append(task)

verbose = {
	"name" : "verbose",
	"doc" : {
		"brief" : "Verbose modes",
	},
	"values" : {
		0 : {
    		"name": "PastixVerboseNot",
    		"value": "0",
    		"brief": "Nothing"
    	},
		1 : {
    		"name": "PastixVerboseNo",
    		"value": "1",
    		"brief": "Default"
    	},
		2 : {
    		"name": "PastixVerboseYes",
    		"value": "2",
    		"brief": "Extended"
    	}
	}
}

enums.append(verbose)

io = {
	"name" : "io",
	"doc" : {
		"brief" : "IO strategy for graph and ordering",
	},
	"values" : {
		0 : {
    		"name": "PastixIONo",
    		"value": "0",
    		"brief": "No output or input"
    	},
		1 : {
    		"name": "PastixIOLoad",
    		"value": "1",
    		"brief": "Load ordering and symbol matrix instead of applying symbolic factorization step"
    	},
		2 : {
    		"name": "PastixIOSave",
    		"value": "2",
    		"brief": "Save ordering and symbol matrix after symbolic factorization step"
    	},
		3 : {
    		"name": "PastixIOLoadGraph",
    		"value": "4",
    		"brief": "Load graph  during ordering step"
    	},
		4 : {
    		"name": "PastixIOSaveGraph",
    		"value": "8",
    		"brief": "Save graph  during ordering step"
    	},
		5 : {
    		"name": "PastixIOLoadCSC",
    		"value": "16",
    		"brief": "Load CSC(d) during ordering step"
    	},
		6 : {
    		"name": "PastixIOSaveCSC",
    		"value": "32",
    		"brief": "Save CSC(d) during ordering step"
    	}
	}
}

enums.append(io)

fact_mode = {
	"name" : "fact_mode",
	"doc" : {
		"brief" : "Factorization Schur modes",
		"details" : r'''
Describe which part of the matrix is factorized or not
''',
	},
	"values" : {
		0 : {
    		"name": "PastixFactModeLocal",
    		"value": "0"
    	},
		1 : {
    		"name": "PastixFactModeSchur",
    		"value": "1"
    	},
		2 : {
    		"name": "PastixFactModeBoth",
    		"value": "2"
    	}
	}
}

enums.append(fact_mode)

solv_mode = {
	"name" : "solv_mode",
	"doc" : {
		"brief" : "Solve Schur modes",
		"details" : r'''
Describe which part of the solve is applied with the matrix

\f[ A = \left( \begin{array}{cc}
            L_{11}U_{11} & U_{12} \\
            L_{21}       & S_{22} \end{array} \right) \f]

For the lower part (and symmetrically for upper part):
  -# Solve \f[ L_{11} * x_{11} = b_{11} \f]
  -# Apply the update \f[ b_{22} = b_{22} - L_{21} * b_{11} \f]
  -# Solve the lower part of \f[ S_{22} * x_{22} = b_{22} \f] if S22 has been previously factorized.

PastixSolvModeLocal applies only the step 1.
PastixSolvModeInterface applies steps 1 and 2.
PastixSolvModeSchur applies all steps.
''',
	},
	"values" : {
		0 : {
    		"name": "PastixSolvModeLocal",
    		"value": "0"
    	},
		1 : {
    		"name": "PastixSolvModeInterface",
    		"value": "1"
    	},
		2 : {
    		"name": "PastixSolvModeSchur",
    		"value": "2"
    	}
	}
}

enums.append(solv_mode)

refine = {
	"name" : "refine",
	"doc" : {
		"brief" : "Iterative refinement algorithms",
	},
	"values" : {
		0 : {
    		"name": "PastixRefineGMRES",
    		"brief": "GMRES"
    	},
		1 : {
    		"name": "PastixRefineCG",
    		"brief": "Conjugate Gradient"
    	},
		2 : {
    		"name": "PastixRefineSR",
    		"brief": "Simple refinement"
    	},
		3 : {
    		"name": "PastixRefineBiCGSTAB",
    		"brief": "BiCGStab"
    	}
	}
}

enums.append(refine)

coeftype = {
	"name" : "coeftype",
	"doc" : {
		"brief" : "Arithmetic types.",
		"details" : r'''
This describes the different arithmetics that can be stored in a sparse matrix.
@remark The values start at 2 for compatibility purpose with PLASMA and
DPLASMA libraries, and they match the ones used in spm.

@sa spm_coeftype_t

@{
''',
	},
	"values" : {
		0 : {
    		"name": "PastixPattern",
    		"value": "SpmPattern"
    	},
		1 : {
    		"name": "PastixFloat",
    		"value": "SpmFloat"
    	},
		2 : {
    		"name": "PastixDouble",
    		"value": "SpmDouble"
    	},
		3 : {
    		"name": "PastixComplex32",
    		"value": "SpmComplex32"
    	},
		4 : {
    		"name": "PastixComplex64",
    		"value": "SpmComplex64"
    	}
	}
}

enums.append(coeftype)

factotype = {
	"name" : "factotype",
	"doc" : {
		"brief" : "Factorization algorithms available for IPARM_FACTORIZATION parameter",
	},
	"values" : {
		0 : {
    		"name": "PastixFactPOTRF",
    		"value": "0",
    		"brief": "Cholesky factorization"
    	},
		1 : {
    		"name": "PastixFactSYTRF",
    		"value": "1",
    		"brief": "LDL^t factorization"
    	},
		2 : {
    		"name": "PastixFactGETRF",
    		"value": "2",
    		"brief": "LU factorization"
    	},
		3 : {
    		"name": "PastixFactPXTRF",
    		"value": "3",
    		"brief": "LL^t factorization for complex matrices"
    	},
		4 : {
    		"name": "PastixFactHETRF",
    		"value": "4",
    		"brief": "LDL^h factorization for complex matrices"
    	},
		5 : {
    		"name": "PastixFactLLH",
    		"value": "0",
    		"brief": "LL^h factorization for complex matrices"
    	},
		6 : {
    		"name": "PastixFactLDLT",
    		"value": "1",
    		"brief": "LDL^t factorization"
    	},
		7 : {
    		"name": "PastixFactLU",
    		"value": "2",
    		"brief": "LU factorization"
    	},
		8 : {
    		"name": "PastixFactLLT",
    		"value": "3",
    		"brief": "LL^t factorization"
    	},
		9 : {
    		"name": "PastixFactLDLH",
    		"value": "4",
    		"brief": "LDL^h factorization for complex matrices"
    	}
	}
}

enums.append(factotype)

scheduler = {
	"name" : "scheduler",
	"doc" : {
		"brief" : "Scheduler",
	},
	"values" : {
		0 : {
    		"name": "PastixSchedSequential",
    		"value": "0",
    		"brief": "Sequential"
    	},
		1 : {
    		"name": "PastixSchedStatic",
    		"value": "1",
    		"brief": "Shared memory with static scheduler"
    	},
		2 : {
    		"name": "PastixSchedParsec",
    		"value": "2",
    		"brief": "PaRSEC scheduler"
    	},
		3 : {
    		"name": "PastixSchedStarPU",
    		"value": "3",
    		"brief": "StarPU scheduler"
    	},
		4 : {
    		"name": "PastixSchedDynamic",
    		"value": "4",
    		"brief": "Shared memory with dynamic scheduler"
    	}
	}
}

enums.append(scheduler)

ordering = {
	"name" : "ordering",
	"doc" : {
		"brief" : "Ordering strategy",
	},
	"values" : {
		0 : {
    		"name": "PastixOrderScotch",
    		"brief": "Use Scotch ordering"
    	},
		1 : {
    		"name": "PastixOrderMetis",
    		"brief": "Use Metis ordering"
    	},
		2 : {
    		"name": "PastixOrderPersonal",
    		"brief": "Apply user's permutation, or load from file"
    	},
		3 : {
    		"name": "PastixOrderPtScotch",
    		"brief": "Use Pt-Scotch ordering"
    	},
		4 : {
    		"name": "PastixOrderParMetis",
    		"brief": "Use ParMetis ordering"
    	}
	}
}

enums.append(ordering)

mpithreadmode = {
	"name" : "mpithreadmode",
	"doc" : {
		"brief" : "MPI thread mode",
	},
	"values" : {
		0 : {
    		"name": "PastixMpiNone",
    		"value": "0",
    		"brief": "No MPI support"
    	},
		1 : {
    		"name": "PastixMpiThreadSingle",
    		"value": "1",
    		"brief": "MPI thread single support"
    	},
		2 : {
    		"name": "PastixMpiThreadFunneled",
    		"value": "2",
    		"brief": "MPI thread funneled support"
    	},
		3 : {
    		"name": "PastixMpiThreadSerialized",
    		"value": "3",
    		"brief": "MPI thread serialized support"
    	},
		4 : {
    		"name": "PastixMpiThreadMultiple",
    		"value": "4",
    		"brief": "MPI thread multiple support"
    	}
	}
}

enums.append(mpithreadmode)

error = {
	"name" : "error",
	"doc" : {
		"brief" : "Error codes",
	},
	"values" : {
		0 : {
    		"name": "PASTIX_SUCCESS",
    		"value": "0",
    		"brief": "No error"
    	},
		1 : {
    		"name": "PASTIX_ERR_UNKNOWN",
    		"value": "1",
    		"brief": "Unknown error"
    	},
		2 : {
    		"name": "PASTIX_ERR_ALLOC",
    		"value": "2",
    		"brief": "Allocation error"
    	},
		3 : {
    		"name": "PASTIX_ERR_NOTIMPLEMENTED",
    		"value": "3",
    		"brief": "Not implemented feature"
    	},
		4 : {
    		"name": "PASTIX_ERR_OUTOFMEMORY",
    		"value": "4",
    		"brief": "Not enough memory"
    	},
		5 : {
    		"name": "PASTIX_ERR_THREAD",
    		"value": "5",
    		"brief": "Error with threads"
    	},
		6 : {
    		"name": "PASTIX_ERR_INTERNAL",
    		"value": "6",
    		"brief": "Internal error"
    	},
		7 : {
    		"name": "PASTIX_ERR_BADPARAMETER",
    		"value": "7",
    		"brief": "Bad parameters given"
    	},
		8 : {
    		"name": "PASTIX_ERR_FILE",
    		"value": "8",
    		"brief": "Error in In/Out operations"
    	},
		9 : {
    		"name": "PASTIX_ERR_INTEGER_TYPE",
    		"value": "9",
    		"brief": "Error with integer types"
    	},
		10 : {
    		"name": "PASTIX_ERR_IO",
    		"value": "10",
    		"brief": "Error with input/output"
    	},
		11 : {
    		"name": "PASTIX_ERR_MPI",
    		"value": "11",
    		"brief": "Error with MPI calls"
    	}
	}
}

enums.append(error)

compress_when = {
	"name" : "compress_when",
	"doc" : {
		"brief" : "Compression strategy available for IPARM_COMPRESS_WHEN parameter",
	},
	"values" : {
		0 : {
    		"name": "PastixCompressNever",
    		"brief": "Do not use compression"
    	},
		1 : {
    		"name": "PastixCompressWhenBegin",
    		"brief": "Compress before any numerical operation (Minimal-Memory)"
    	},
		2 : {
    		"name": "PastixCompressWhenEnd",
    		"brief": "Compress after contributions were accumulated (Just-In-Time)"
    	},
		3 : {
    		"name": "PastixCompressWhenDuring",
    		"brief": "Compress after contributions from other supernodes were accumulated"
    	}
	}
}

enums.append(compress_when)

compress_method = {
	"name" : "compress_method",
	"doc" : {
		"brief" : "Compression method available for IPARM_COMPRESS_METHOD parameter",
	},
	"values" : {
		0 : {
    		"name": "PastixCompressMethodSVD",
    		"brief": "Use singular value decomposition for low-rank compression"
    	},
		1 : {
    		"name": "PastixCompressMethodPQRCP",
    		"brief": "Use partial QR with column pivoting for low-rank compression"
    	},
		2 : {
    		"name": "PastixCompressMethodRQRCP",
    		"brief": "Use randomized QR with column pivoting for low-rank compression"
    	},
		3 : {
    		"name": "PastixCompressMethodTQRCP",
    		"brief": "Use truncated QR with column pivotingfor low-rank compression"
    	},
		4 : {
    		"name": "PastixCompressMethodRQRRT",
    		"brief": "Use randomized QR with rotation for low-rank compression"
    	},
		5 : {
    		"name": "PastixCompressMethodNbr",
    		"brief": "Total number of available compression methods"
    	}
	}
}

enums.append(compress_method)

compress_ortho = {
	"name" : "compress_ortho",
	"doc" : {
		"brief" : "Orthogonalization method available for IPARM_COMPRESS_ORTHO parameter",
	},
	"values" : {
		0 : {
    		"name": "PastixCompressOrthoCGS",
    		"brief": "Orthogonalize low-rank bases with Gram-Schimdt"
    	},
		1 : {
    		"name": "PastixCompressOrthoQR",
    		"brief": "Orthogonalize low-rank bases with QR decomposition"
    	},
		2 : {
    		"name": "PastixCompressOrthoPartialQR",
    		"brief": "Orthogonalize low-rank bases with projections in orthogonal space followed by smaller QR"
    	}
	}
}

enums.append(compress_ortho)

split = {
	"name" : "split",
	"doc" : {
		"brief" : "Splitting strategy available for IPARM_SPLITTING_STRATEGY parameter",
	},
	"values" : {
		0 : {
    		"name": "PastixSplitNot",
    		"brief": "Do not apply dedicated low-rank clustering strategy"
    	},
		1 : {
    		"name": "PastixSplitKway",
    		"brief": "Use k-way partitioning"
    	},
		2 : {
    		"name": "PastixSplitKwayProjections",
    		"brief": "Use projections and k-way in clusters"
    	}
	}
}

enums.append(split)

layout = {
	"name" : "layout",
	"doc" : {
		"brief" : "Direction of the matrix storage",
	},
	"values" : {
		0 : {
    		"name": "PastixRowMajor",
    		"value": "101",
    		"brief": "Storage in row major order"
    	},
		1 : {
    		"name": "PastixColMajor",
    		"value": "102",
    		"brief": "Storage in column major order"
    	}
	}
}

enums.append(layout)

trans = {
	"name" : "trans",
	"doc" : {
		"brief" : "Transpostion",
	},
	"values" : {
		0 : {
    		"name": "PastixNoTrans",
    		"value": "111",
    		"brief": "Use A"
    	},
		1 : {
    		"name": "PastixTrans",
    		"value": "112",
    		"brief": "Use A^t"
    	},
		2 : {
    		"name": "PastixConjTrans",
    		"value": "113",
    		"brief": "Use conj(A^t)"
    	}
	}
}

enums.append(trans)

mtxtype = {
	"name" : "mtxtype",
	"doc" : {
		"brief" : "Matrix symmetry type property.",
		"details" : r'''
@remark Must match transposition.
''',
	},
	"values" : {
		0 : {
    		"name": "PastixGeneral",
    		"value": "PastixNoTrans",
    		"brief": "The matrix is general"
    	},
		1 : {
    		"name": "PastixSymmetric",
    		"value": "PastixTrans",
    		"brief": "The matrix is symmetric"
    	},
		2 : {
    		"name": "PastixHermitian",
    		"value": "PastixConjTrans",
    		"brief": "The matrix is hermitian"
    	}
	}
}

enums.append(mtxtype)

uplo = {
	"name" : "uplo",
	"doc" : {
		"brief" : "Upper/Lower part",
	},
	"values" : {
		0 : {
    		"name": "PastixUpper",
    		"value": "121",
    		"brief": "Use lower triangle of A"
    	},
		1 : {
    		"name": "PastixLower",
    		"value": "122",
    		"brief": "Use upper triangle of A"
    	},
		2 : {
    		"name": "PastixUpperLower",
    		"value": "123",
    		"brief": "Use the full A"
    	}
	}
}

enums.append(uplo)

coefside = {
	"name" : "coefside",
	"doc" : {
		"brief" : "Data blocks used in the kernel",
		"details" : r'''
@warning Must be 0 and 1 respectively for Left and Upper as is it used to
shift the pointers in the kernels from the lower to upper part.
''',
	},
	"values" : {
		0 : {
    		"name": "PastixLCoef",
    		"value": "0",
    		"brief": "Coefficients of the lower triangular L are used"
    	},
		1 : {
    		"name": "PastixUCoef",
    		"value": "1",
    		"brief": "Coefficients of the upper triangular U are used"
    	},
		2 : {
    		"name": "PastixLUCoef",
    		"value": "2",
    		"brief": "Coefficients of the upper/lower triangular U/L are used"
    	}
	}
}

enums.append(coefside)

diag = {
	"name" : "diag",
	"doc" : {
		"brief" : "Diagonal",
	},
	"values" : {
		0 : {
    		"name": "PastixNonUnit",
    		"value": "131",
    		"brief": "Diagonal is non unitary"
    	},
		1 : {
    		"name": "PastixUnit",
    		"value": "132",
    		"brief": "Diagonal is unitary"
    	}
	}
}

enums.append(diag)

side = {
	"name" : "side",
	"doc" : {
		"brief" : "Side of the operation",
	},
	"values" : {
		0 : {
    		"name": "PastixLeft",
    		"value": "141",
    		"brief": "Apply operator on the left"
    	},
		1 : {
    		"name": "PastixRight",
    		"value": "142",
    		"brief": "Apply operator on the right"
    	}
	}
}

enums.append(side)

normtype = {
	"name" : "normtype",
	"doc" : {
		"brief" : "Norms",
	},
	"values" : {
		0 : {
    		"name": "PastixOneNorm",
    		"value": "171",
    		"brief": "One norm:       max_j( sum_i( |a_{ij}| ) )"
    	},
		1 : {
    		"name": "PastixFrobeniusNorm",
    		"value": "174",
    		"brief": "Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) )"
    	},
		2 : {
    		"name": "PastixInfNorm",
    		"value": "175",
    		"brief": "Inifinite norm: max_i( sum_j( |a_{ij}| ) )"
    	},
		3 : {
    		"name": "PastixMaxNorm",
    		"value": "177",
    		"brief": "Inifinite norm: max_{i,j}( | a_{ij} | )"
    	}
	}
}

enums.append(normtype)

dir = {
	"name" : "dir",
	"doc" : {
		"brief" : "Direction",
	},
	"values" : {
		0 : {
    		"name": "PastixDirForward",
    		"value": "391",
    		"brief": "Forward direction"
    	},
		1 : {
    		"name": "PastixDirBackward",
    		"value": "392",
    		"brief": "Backward direction"
    	}
	}
}

enums.append(dir)

