#=

 @file PaStiX.jl

 PaStiX julia wrapper

 @copyright 2020-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.2.1
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date 2022-10-11

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_julia

=#

module PaStiX
using CBinding
using Libdl
include("pastix_enums.jl")

function pastix_library_path()
    x = Libdl.dlext
    return "libpastix.$x"
end
libpastix = pastix_library_path()

@cbindings libpastix begin
    @cextern pastixOrderInit( ordeptr::Ptr{Pastix_order_t}, baseval::Pastix_int_t, vertnbr::Pastix_int_t, cblknbr::Pastix_int_t, perm::Ptr{Pastix_int_t}, invp::Ptr{Pastix_int_t}, rang::Ptr{Pastix_int_t}, tree::Ptr{Pastix_int_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderAlloc( ordeptr::Ptr{Pastix_order_t}, vertnbr::Pastix_int_t, cblknbr::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderAllocId( ordeptr::Ptr{Pastix_order_t}, vertnbr::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderExit( ordeptr::Ptr{Pastix_order_t} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixOrderBase( ordeptr::Ptr{Pastix_order_t}, baseval::Pastix_int_t )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixOrderCheck( ordeptr::Ptr{Pastix_order_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderCopy( ordedst::Ptr{Pastix_order_t}, ordesrc::Ptr{Pastix_order_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderGet( pastix_data::Ptr{Pastix_data_t} )::Ptr{Pastix_order_t}
end

@cbindings libpastix begin
    @cextern pastixOrderBcast( ordemesh::Ptr{Pastix_order_t}, root::Cint, pastix_comm::__get_mpi_type__() )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixOrderGrid( myorder::Ptr{Cvoid}, nx::Pastix_int_t, ny::Pastix_int_t, nz::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderLoad( pastix_data::Ptr{Pastix_data_t}, ordeptr::Ptr{Pastix_order_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastixOrderSave( pastix_data::Ptr{Pastix_data_t}, ordeptr::Ptr{Pastix_order_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix( pastix_data::Ptr{Cvoid}, pastix_comm::__get_mpi_type__(), n::Pastix_int_t, colptr::Ptr{Pastix_int_t}, rowptr::Ptr{Pastix_int_t}, values::Ptr{Cvoid}, perm::Ptr{Pastix_int_t}, invp::Ptr{Pastix_int_t}, B::Ptr{Cvoid}, nrhs::Pastix_int_t, iparm::Ptr{Pastix_int_t}, dparm::Ptr{Cdouble} )::Cint
end

@cbindings libpastix begin
    @cextern pastixInitParam( iparm::Ptr{Pastix_int_t}, dparm::Ptr{Cdouble} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixInit( pastix_data::Ptr{Cvoid}, pastix_comm::__get_mpi_type__(), iparm::Ptr{Pastix_int_t}, dparm::Ptr{Cdouble} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixInitWithAffinity( pastix_data::Ptr{Cvoid}, pastix_comm::__get_mpi_type__(), iparm::Ptr{Pastix_int_t}, dparm::Ptr{Cdouble}, bindtab::Ptr{Cint} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixFinalize( pastix_data::Ptr{Cvoid} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastix_task_analyze( pastix_data::Ptr{Pastix_data_t}, spm::Ptr{spm.spmatrix_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_task_numfact( pastix_data::Ptr{Pastix_data_t}, spm::Ptr{spm.spmatrix_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_task_solve( pastix_data::Ptr{Pastix_data_t}, m::Pastix_int_t, nrhs::Pastix_int_t, B::Ptr{Cvoid}, ldb::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_task_refine( pastix_data::Ptr{Pastix_data_t}, n::Pastix_int_t, nrhs::Pastix_int_t, B::Ptr{Cvoid}, ldb::Pastix_int_t, X::Ptr{Cvoid}, ldx::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_order( pastix_data::Ptr{Pastix_data_t}, spm::Ptr{spm.spmatrix_t}, myorder::Ptr{Pastix_order_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_symbfact( pastix_data::Ptr{Pastix_data_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_reordering( pastix_data::Ptr{Pastix_data_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_blend( pastix_data::Ptr{Pastix_data_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_spm2bcsc( pastix_data::Ptr{Pastix_data_t}, spm::Ptr{spm.spmatrix_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_bcsc2ctab( pastix_data::Ptr{Pastix_data_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_sopalin( pastix_data::Ptr{Pastix_data_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastixRhsInit( rhs::Ptr{Pastix_rhs_t} )::Cint
end

@cbindings libpastix begin
    @cextern pastixRhsFinalize( rhs::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixRhsDoubletoSingle( dB::Pastix_rhs_t, sB::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixRhsSingleToDouble( sB::Pastix_rhs_t, dB::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_applyorder( pastix_data::Ptr{Pastix_data_t}, dir::spm.spm_dir_t, m::Pastix_int_t, n::Pastix_int_t, B::Ptr{Cvoid}, ldb::Pastix_int_t, Bp::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_trsm( pastix_data::Ptr{Pastix_data_t}, side::Cint, uplo::Cint, trans::Pastix_trans_t, diag::Pastix_diag_t, b::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_diag( pastix_data::Ptr{Pastix_data_t}, b::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_solve( pastix_data::Ptr{Pastix_data_t}, b::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_refine( pastix_data::Ptr{Pastix_data_t}, b::Pastix_rhs_t, x::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastix_subtask_solve_adv( pastix_data::Ptr{Pastix_data_t}, transA::Pastix_trans_t, b::Pastix_rhs_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixSetSchurUnknownList( pastix_data::Ptr{Pastix_data_t}, n::Pastix_int_t, list::Ptr{Pastix_int_t} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixGetSchur( pastix_data::Ptr{Pastix_data_t}, S::Ptr{Cvoid}, lds::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixExpand( pastix_data::Ptr{Pastix_data_t}, spm::Ptr{spm.spmatrix_t} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixGetDiag( pastix_data::Ptr{Pastix_data_t}, x::Ptr{Cvoid}, incx::Pastix_int_t )::Cint
end

@cbindings libpastix begin
    @cextern pastixGetOptions( argc::Cint, argv::Ptr{Cvoid}, iparm::Ptr{Pastix_int_t}, dparm::Ptr{Cdouble}, check::Ptr{Cint}, driver::Ptr{spm.spm_driver_t}, filename::Ptr{Cvoid} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixDumpParam( pastix_data::Ptr{Pastix_data_t} )::Cvoid
end

@cbindings libpastix begin
    @cextern pastixCheckParam( iparm::Ptr{Pastix_int_t}, dparm::Ptr{Cdouble} )::Cint
end

end #module
