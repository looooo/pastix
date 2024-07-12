/**
 *
 * @file parse_options.h
 *
 * This file is generated automatically. If you want to modify it, modify
 * ${PASTIX_HOME}/tools/gen_param/pastix_[iparm/dparm/enums].py and run
 * ${PASTIX_HOME}/tools/gen_param/gen_parm_files.py ${PASTIX_HOME}.
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @author Alycia Lisito
 * @author Tom Moenne-Loccoz
 * @date 2024-07-05
 *
 */
#ifndef _parse_options_h_
#define _parse_options_h_

BEGIN_C_DECLS

pastix_iparm_t parse_iparm( const char *iparm  );
pastix_dparm_t parse_dparm( const char *dparm  );
int            parse_enums( const char *string );

const char* pastix_task_getstr( pastix_task_t value );
const char* pastix_trace_getstr( pastix_trace_t value );
const char* pastix_verbose_getstr( pastix_verbose_t value );
const char* pastix_io_getstr( pastix_io_t value );
const char* pastix_fact_mode_getstr( pastix_fact_mode_t value );
const char* pastix_solv_mode_getstr( pastix_solv_mode_t value );
const char* pastix_refine_getstr( pastix_refine_t value );
const char* pastix_coeftype_getstr( pastix_coeftype_t value );
const char* pastix_factotype_getstr( pastix_factotype_t value );
const char* pastix_factolookside_getstr( pastix_factolookside_t value );
const char* pastix_scheduler_getstr( pastix_scheduler_t value );
const char* pastix_ordering_getstr( pastix_ordering_t value );
const char* pastix_mpithreadmode_getstr( pastix_mpithreadmode_t value );
const char* pastix_error_getstr( pastix_error_t value );
const char* pastix_compress_when_getstr( pastix_compress_when_t value );
const char* pastix_compress_method_getstr( pastix_compress_method_t value );
const char* pastix_compress_ortho_getstr( pastix_compress_ortho_t value );
const char* pastix_split_getstr( pastix_split_t value );
const char* pastix_layout_getstr( pastix_layout_t value );
const char* pastix_trans_getstr( pastix_trans_t value );
const char* pastix_mtxtype_getstr( pastix_mtxtype_t value );
const char* pastix_uplo_getstr( pastix_uplo_t value );
const char* pastix_coefside_getstr( pastix_coefside_t value );
const char* pastix_diag_getstr( pastix_diag_t value );
const char* pastix_side_getstr( pastix_side_t value );
const char* pastix_normtype_getstr( pastix_normtype_t value );
const char* pastix_dir_getstr( pastix_dir_t value );

void pastix_param2csv( const pastix_data_t *pastix_data, FILE *csv );
int iparm_check_values( const pastix_int_t *iparm );
int dparm_check_values( const double       *dparm );

END_C_DECLS

#endif /* _parse_options_h_ */
