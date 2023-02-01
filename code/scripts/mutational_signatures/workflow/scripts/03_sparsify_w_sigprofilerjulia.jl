"""
@created: 22 Jun 2021
@modified: 22 Sep 2021
@author: Yoann Pradat

Run the sparsity-enforcement procedure of sigprofilerjulia algorithm.
"""

# Libraries ============================================================================================================

using Distributed
using Logging

@everywhere using ArgParse
@everywhere using DelimitedFiles
@everywhere using CSV
@everywhere using DataFrames
@everywhere using ProgressMeter

# Functions ============================================================================================================

struct FlushingIO{T}
    io::T
end

function Base.println(io::FlushingIO, args...)
    println(io.io, args...)
    flush(io.io)
end

"""
    convert_array_to_df(basis_name)

Convert array to dataframe.
"""
function convert_array_to_df(arr; convert_type=nothing, row_names=false, col_names=false)
    if row_names
        if col_names
            rows = arr[2:end,1]
            cols = arr[1, 2:end]
            vals = arr[2:end, 2:end]
        else
            rows = arr[1:end,1]
            cols = nothing
            vals = arr[1:end, 2:end]
        end
    else
        row = nothing
        if col_names
            cols = arr[1, 1:end]
            vals = arr[2:end, 1:end]
        else
            cols = nothing
            vals = arr
        end
    end

    if !isnothing(convert_type)
        vals = convert(Matrix{convert_type}, vals)
    end

    if col_names
        df = DataFrame(vals, Symbol.(cols))
    else
        df = DataFrame(vals, :auto)
    end

    if row_names
        insertcols!(df, 1, (:Rows => rows));
    end

    df
end

"""
    convert_df_to_array(basis_name)

Convert dataframe to array. The column are lost as well as the row names encoded in the :Rows column if present.
"""
function convert_df_to_array(df)
    df_copy = deepcopy(df)
    if "Rows" in names(df)
        select!(df_copy, Not(:Rows));
    end
    Matrix(df_copy)
end


"""
    load_signature_profiles(cohort, method, basis_name)

Load the matrix of signature counts for the specified cohort.
"""
function load_signature_profiles(basis_name)
    filepath = joinpath("resources/signatures_cosmic", "signatures_" * basis_name * ".tsv");
    df_W = read_with_rows(filepath, true, delim="\t");
    rename!(df_W, :row_name => :Rows)
end


@everywhere function remove_sigs(j, parallelExposure, parallelProcesses, originalGenomes)
    """
        remove_sig(j)

    Make signatures sparse for j^th sample.
    """
    exposuresOutput, accr, frob_rel_div, norm_one_dif = removeAllSingleSignatures(parallelExposure[:,j],
                                                                                  parallelProcesses,
                                                                                  originalGenomes[:,j])
    return j, exposuresOutput, accr, frob_rel_div, norm_one_dif
end

function main(args)

    @eval @everywhere path_root=args["sigprofilerjulia"]
    include(joinpath(pwd(),path_root,"src/removeAllSingleSignatures.jl"))
    include(joinpath(pwd(),path_root,"src/removeOneSig.jl"))
    include(joinpath(pwd(),path_root,"src/evalSingleSample.jl"))

    df_M = convert_array_to_df(readdlm(args["count_mut"], '\t'), convert_type=Float64, row_names=true, col_names=true);
    df_W = convert_array_to_df(readdlm(args["basis"], '\t'), convert_type=Float64, row_names=true, col_names=true);
    df_H = convert_array_to_df(readdlm(args["count_sig"], '\t'), convert_type=Float64, row_names=true, col_names=true);

    genomes = deepcopy(convert_df_to_array(df_M));
    @everywhere originalGenomes=$genomes

    # signature fitting
    exposure = deepcopy(convert_df_to_array(df_H));
    processes = deepcopy(convert_df_to_array(df_W));

    @eval @everywhere parallelExposure=$exposure
    @eval @everywhere parallelProcesses=$processes

    nSignatures = size(parallelExposure, 1);
    nGenomes = size(parallelExposure, 2);

    # run in parallel
    pmap_res = @showprogress pmap(1:size(parallelExposure, 2)) do j
        remove_sigs(j, parallelExposure, parallelProcesses, originalGenomes)
    end

    exposureSparse = (x -> x[2]).(pmap_res)
    exposureSparse = reshape(vcat(exposureSparse'...), (nGenomes, nSignatures))'

    exposureSparse = hcat(df_H[!,"Rows"], exposureSparse)
    exposureSparse = vcat(reshape(names(df_M), (1, nGenomes+1)), exposureSparse)
    exposureSparse[1,1] = readdlm(args["count_sig"], '\t')[1,1]

    # save output
    writedlm(args["output"], exposureSparse, '\t')
end

# run ==================================================================================================================

@eval @everywhere s = ArgParseSettings()
@everywhere @add_arg_table s begin
    "--count_mut"
        help = "Path to the original matrix of mutation counts."
        default = "../../results/mutational_signatures/counts_mutations/counts_mutations_sbs_96_min_mut_met500.tsv"
        arg_type = String
    "--count_sig"
        help = "Path to the matrix of signature counts."
        default = "../../results/mutational_signatures/projection_known_signatures/MutationalPatterns/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_met500.tsv"
        arg_type = String
    "--basis"
        help = "Path to the matrix of signature profiles."
        default = "resources/signatures_cosmic/signatures_cosmic_sbs_96_v3.2.tsv"
        arg_type = String
    "--sigprofilerjulia"
        help = "Path to sigprofilerjulia root folder"
        default = "external/sigprofilerjulia"
        arg_type = String
    "--output"
        help = "Path to the matrix of signature counts sparsified."
        arg_type = String
end

@eval @everywhere args = parse_args(s);
out = FlushingIO(stderr)

println(out, "Parsed args:")
for (arg,val) in args
    println(out, "  $arg  =>  $val")
end

main(args)
