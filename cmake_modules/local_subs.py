"""
 @file local_subs.py

 Python PaStiX specific substitution rules for the Precision Generator script.

 @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.4.0
 @author Mathieu Faverge
 @author Esragul Korkmaz
 @author Tony Delarue
 @author Alycia Lisito
 @author Brieuc Nicolas
 @date 2023-12-18

"""
subs = {
    # ------------------------------------------------------------
    # replacements applied to normal precision files.
    'normal': [
        # pattern                single                  double                  single-complex          double-complex
        #'12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890')

        # ----- Additional Datatypes
        ('int',                  'float',                'double',               'spm_complex32_t',     r'\bspm_complex64_t'   ),
        ('SpmPattern',           'SpmFloat',             'SpmDouble',            'SpmComplex32',        r'\bSpmComplex64'      ),
        ('SpmPattern',           'SpmFloat',             'SpmDouble',            'SpmFloat',            r'\bSpmDouble'         ),
        ('int',                  'float',                'double',               'pastix_complex32_t',  r'\bpastix_complex64_t'),
        ('int',                  'float',                'float',                'pastix_complex32_t',  r'\bpastix_complex32_t'),

        ('PastixPattern',        'PastixFloat',          'PastixDouble',         'PastixComplex32',     r'\bPastixComplex64'   ),
        ('PastixPattern',        'PastixFloat',          'PastixFloat',          'PastixComplex32',     r'\bPastixComplex32'   ),
        ('PastixPattern',        'PastixFloat',          'PastixDouble',         'PastixFloat',         r'\bPastixDouble'      ),
        ('MPI_INT',              'MPI_FLOAT',            'MPI_DOUBLE',           'MPI_COMPLEX32',        'MPI_COMPLEX64'       ),

        # ----- Additional PaStiX BLAS
        ('',                     '_slr',                 '_dlr',                 '_clr',                 '_zlr'                ),
        ('',                     'sgelrops',             'dgelrops',             'cgelrops',             'zgelrops'            ),
        ('',                     'sge2lr',               'dge2lr',               'cge2lr',               'zge2lr'              ),
        ('',                     'srradd',               'drradd',               'crradd',               'zrradd'              ),
        ('',                     'spqrcp',               'dpqrcp',               'cpqrcp',               'zpqrcp'              ),
        ('',                     'sdiag' ,               'ddiag' ,               'cdiag' ,               'zdiag'               ),
        ('',                     'spotrf',               'dpotrf',               'cpxtrf',               'zpxtrf'              ),
        ('',                     'srqrcp',               'drqrcp',               'crqrcp',               'zrqrcp'              ),
        ('',                     'srqrrt',               'drqrrt',               'crqrrt',               'zrqrrt'              ),
        ('',                     'stqrcp',               'dtqrcp',               'ctqrcp',               'ztqrcp'              ),
        ('',                     'sxx2fr',               'dxx2fr',               'cxx2fr',               'zxx2fr'              ),
        ('',                     'sxx2lr',               'dxx2lr',               'cxx2lr',               'zxx2lr'              ),
        ('',                     'sytrf',                'sytrf',                'hetrf',                'hetrf'               ),
        ('',                     'slassq',               'dlassq',               'slassq',               'dlassq'              ),

        # ----- PaStiX Variables
        (r'\b',                 r'szero\b',             r'dzero\b',             r'czero\b',             r'zzero\b'             ),
        (r'\b',                 r'sone\b',              r'done\b',              r'cone\b',              r'zone\b'              ),

        # ----- SPM Prefixes
        ('spm_p',                'spm_s',                'spm_d',                'spm_c',                'spm_z'               ),

        # ----- PaStiX Prefixes
        ('CORE_P',               'CORE_S',               'CORE_D',               'CORE_C',               'CORE_Z'              ),
        ('blok_p',               'blok_s',               'blok_d',               'blok_c',               'blok_z'              ),
        ('cblk_p',               'cblk_s',               'cblk_d',               'cblk_c',               'cblk_z'              ),
        ('coeftab_p',            'coeftab_s',            'coeftab_d',            'coeftab_c',            'coeftab_z'           ),
        ('cuda_p',               'cuda_s',               'cuda_d',               'cuda_c',               'cuda_z'              ),
        ('core_p',               'core_s',               'core_d',               'core_c',               'core_z'              ),
        ('csc_p',                'csc_s',                'csc_d',                'csc_c',                'csc_z'               ),
        ('pastix_p',             'pastix_s',             'pastix_d',             'pastix_c',             'pastix_z'            ),
        ('sequential_p',         'sequential_s',         'sequential_d',         'sequential_c',         'sequential_z'        ),
        ('thread_p',             'thread_s',             'thread_d',             'thread_c',             'thread_z'            ),
        ('static_p',             'static_s',             'static_d',             'static_c',             'static_z'            ),
        ('dynamic_p',            'dynamic_s',            'dynamic_d',            'dynamic_c',            'dynamic_z'           ),
        ('runtime_p',            'runtime_s',            'runtime_d',            'runtime_c',            'runtime_z'           ),
        ('vec_p',                'vec_s',                'vec_d',                'vec_c',                'vec_z'               ),
        ('task_p',               'task_s',               'task_d',               'task_c',               'task_z'              ),
        ('solve_p',              'solve_s',              'solve_d',              'solve_c',              'solve_z'             ),
        ('',                     'slag2d',               'slag2d',               'clag2z',               'clag2z'              ),
    ], #end normal

    # ------------------------------------------------------------
    # replacements applied to mixed precision files.
    'mixed' : [
        # double/single,          double/single-complex
        #'12345678901234567890', '12345678901234567890')

        # ----- Additional Datatypes
        ('double',               'pastix_complex64_t'  ),
        ('float',                'pastix_complex32_t'  ),
        ('double',               'spm_complex64_t'     ),

        ('coeftab_dcblk',        'coeftab_zcblk'       ),
        ('coeftab_ds',           'coeftab_zc'          ),
        ('cpucblk_ds',           'cpucblk_zc'          ),
        ('pastix_dcores',        'pastix_zcores'       ),
        ('pastix_ds',            'pastix_zc'           ),
        ('pastix_scores',        'pastix_ccores'       ),
        ('cpucblk_s',            'cpucblk_c'           ),
        ('DSOVERFCHECK',         'DOVERFCHECK'         ),
        ('DSOVERFCHECK',         'ZCOVERFCHECK'        ),
    ] #end mixed
}

exceptfrom = []
