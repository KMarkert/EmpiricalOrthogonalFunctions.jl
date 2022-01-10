
"""
Extract the principal component time series (PCs)
"""
function pcs(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
    scaling::Symbol = :none,
)
    if isnothing(n)
        n = length(eof.eigenvals)
    end

    out = eof.pcs[:, 1:n]

    if scaling == :divide
        # Divide by the square-root of the eigenvalue.
        for i = 1:n
            out[:, i] = out[:, i] ./ sqrt(eof.eigenvals[i])
        end
    elseif scaling == :multiply
        # Multiply by the square root of the eigenvalue.
        for i = 1:n
            out[:, i] = out[:, i] .* sqrt(eof.eigenvals[i])
        end
    elseif scaling != :none
        @error(
            "Could not understand scaling option. Select either :none, :divide, or :multiply"
        )
        return
    end

    return out
end

"""
Extract the empirical orthogonal functions (EOFs)
"""
function eofs(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
    scaling::Symbol = :none,
)
    if isnothing(n)
        n = length(eof.eigenvals)
    end

    out = eof.eofs[:, 1:n]

    if scaling == :divide
        # Divide by the square-root of the eigenvalue.
        for i = 1:n
            out[:, i] = out[:, i] ./ sqrt(eof.eigenvals[i])
        end
    elseif scaling == :multiply
        # Multiply by the square root of the eigenvalue.
        for i = 1:n
            out[:, i] = out[:, i] .* sqrt(eof.eigenvals[i])
        end
    elseif scaling != :none
        @error(
            "Could not understand scaling option. Select either :none, :divide, or :multiply"
        )
        return
    end

    return out
end

"""
Empirical orthogonal functions (EOFs) expressed as the
correlation between the principal component time series (PCs)
and the time series of the `eof` input *dataset* at each grid
point.
"""
function correlationmap(
    eof::EmpiricalOrthogonalFunction;
    outshape::Union{Nothing,Tuple{Int,Int}} = nothing,
    n::Union{Nothing,Int} = nothing,
)

    p = pcs(eof, n = n, scaling = :divide)

    p_center = Statistics.mean(p, dims = 2)

    pcs_cent = p .- p_center
    field_cent = eof.dataset

    pcs_std = Statistics.std(pcs_cent, dims = 1)
    field_std = field_std = Statistics.std(field_cent, dims = 1)

    div = size(pcs_cent, 1)

    cor = transpose(field_cent' * pcs_cent) / div

    cor = cor ./ (pcs_std' * field_std)

    if !isnothing(outshape)
        cor = reshape(cor, (n, outshape...))
    end

    return cor
end

"""
Empirical orthogonal functions (EOFs) expressed as the
covariance between the principal component time series (PCs)
and the time series of the `eof` input *dataset* at each grid
point.
"""
function covariancemap(
    eof::EmpiricalOrthogonalFunction;
    outshape::Union{Nothing,Tuple{Int,Int}} = nothing,
    n::Union{Nothing,Int} = nothing,
    ddof::Int = 1,
    scaling = :none,
)
    p = pcs(eof, n = n, scaling = scaling)

    p_center = Statistics.mean(p, dims = 2)

    pcs_cent = p .- p_center
    field_cent = eof.dataset

    div = size(pcs_cent, 1) - ddof

    cor = transpose(field_cent' * pcs_cent) / div

    if !isnothing(outshape)
        cor = reshape(cor, (n, outshape...))
    end

    return cor
end

"""
Extract the eigenvalues (decreasing variances) associated with each EOF.
"""
function eigenvalues(eof::EmpiricalOrthogonalFunction; n::Union{Nothing,Int} = nothing)
    if isnothing(n)
        n = length(eof.eigenvals)
    end
    return eof.eigenvals[1:n]
end

"""
Fractional EOF mode variances.
The fraction of the total variance explained by each EOF mode,
values between 0 and 1 inclusive.
"""
function variancefraction(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
)
    return eigenvalues(eof; n = n) / totalanomalyvar(eof)
end

"""
Total variance associated with the field of anomalies (the sum
of the eigenvalues).
"""
function totalanomalyvar(eof::EmpiricalOrthogonalFunction)
    return sum(eof.eigenvals)
end

"""
The method of North et al. (1982) is used to compute the typical
error for each eigenvalue. It is assumed that the number of
times in the input data set is the same as the number of
independent realizations. If this assumption is not valid then
the result may be inappropriate.
"""
function northtest(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
    scaled::Bool = false,
)

    nrecords = size(eof.dataset, 1)
    factor = sqrt(2.0 / nrecords)

    if scaled
        factor /= totalanomalyvar(eof)
    end

    return eigenvalues(eof; n = n) * factor
end

"""
Project a field onto the EOFs.
Given a data set, projects it onto the EOFs to generate a
corresponding set of pseudo-PCs.
"""
function projectfield(eof_arr::Array{<:Any,2}, field::Array{<:Any,2})
    # Project the data set onto the EOFs using a matrix multiplication.
    return field * eof_arr
end

function projectfield(
    eof::EmpiricalOrthogonalFunction,
    field::Array{<:Any,2};
    n::Union{Nothing,Int} = nothing,
    scaling::Symbol = :none,
)

    # extract out eofs for projection
    eofs_arr = eofs(eof, n = n, scaling = scaling)

    # check to make sure all data with valid index in EOF is valid in the input field
    @assert all(!ismissing(field[:, eof.valididx]))

    return projectfield(eofs_arr, field)
end

function projectfield(
    eof::EmpiricalOrthogonalFunction,
    field::Array{<:Any,3};
    n::Union{Nothing,Int} = nothing,
    scaling::Symbol = :none,
    timedim::Int = 3,
)

    geodims = [x for x in 1:3 if x != timedim]
    shape = size(dataset)
    spatial_shape = shape[geodims]
    time_shape = shape[timedim]

    if timedim == 2
        field = permutedims(field, (2, 1, 3))
    elseif timedim == 3
        field = permutedims(flattened, (3, 1, 2))
    end

    flat_field = reshape(field, (time_shape, prod(spatial_shape)))

    projected_pcs = projectfield(eof, flat_field; n = n, scaling = scaling)

    # reshape spatial field back to the original
    return reshape(projected_pcs, (spatial_shape..., time_shape))

end

"""
Reconstructed input data field based on a subset of EOFs.
"""
function reconstruct(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
    uncenter::Bool = true,
)

    eofs_arr = eofs(eof, n = n)
    ps = pcs(eof, n = n)

    rval = (ps * eofs_arr')

    if uncenter
        rval = rval .+ eof.center
    end

    return permutedims(rval, (2, 1))
end


function orthorotation(
    components::Array{<:Any,2};
    method::Symbol = :varimax,
    maxiter::Int = 50,
    tol::Float64 = 1e-6,
)
    nrows, ncols = size(components)

    rotation_matrix = I(ncols)
    var = 0
    i = 0

    for _ = 1:maxiter
        comp_rot = components * rotation_matrix
        if method == :varimax
            tmp = comp_rot * transpose(sum((comp_rot .^ 2), dims = 1) / nrows)
        elseif method == :quartimax
            tmp = 0
        end
        x = @. comp_rot^3 - tmp
        u, s, v = svd(transpose(components) * x)
        rotation_matrix = u * v
        var_new = sum(s)
        if var != 0 && var_new < var * (1 + tol)
            break
        end
        var = var_new
        i += 1
    end

    if i == maxiter
        @warn(
            "matrix rotation did not converge, try specifying a larger max iteration value for stable result"
        )
    end

    return components * rotation_matrix
end

"""
Apply orthogonal rotation to EOF data.
After rotation the original dataset will be projected on the rotated EOF
to create new PCs. Additionally new EOFs and PCs are ordered in decreasing
variance
"""
function orthorotation(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
    method::Symbol = :varimax,
    maxiter::Int = 50,
    tol::Float64 = 1e-6,
)
    components = eofs(eof, n = n)

    valididx = eof.valididx

    reofs = Array{Union{Missing,Float64}}(missing, size(components)...)
    reofs[valididx, :] = orthorotation(
        components[valididx, :];
        method = method,
        maxiter = maxiter,
        tol = tol,
    )

    projected_pcs = projectfield(reofs[valididx, :], eof.dataset[:, valididx])

    rotatedvar = Statistics.var(projected_pcs, dims = 1)[1, :]
    rotatedvarnorm =
        (rotatedvar ./ sum(rotatedvar)) .* sum(variancefraction(eof, n = n))

    sortedidx = sortperm(rotatedvarnorm)[end:-1:1]
    rotatedvar = rotatedvar[sortedidx]

    sorted_reofs = reofs[:, sortedidx]
    sorted_pcs = projected_pcs[:, sortedidx]

    return EmpiricalOrthogonalFunction(
        eof.dataset,
        eof.center,
        eof.valididx,
        sorted_reofs,
        rotatedvar,
        sorted_pcs,
    )
end

"""
Inplace orthorotation that will override the data structure of the input EOF object
"""
function orthorotation!(
    eof::EmpiricalOrthogonalFunction;
    n::Union{Nothing,Int} = nothing,
    method::Symbol = :varimax,
    maxiter::Int = 50,
    tol::Float64 = 1e-6,
)
    if isnothing(n)
        n = length(eof.eigenvals)
    else
        truncate!(eof, n)
    end

    valididx = eof.valididx

    eof.eofs[valididx, :] = orthorotation(
        eof.eofs[valididx, :];
        method = method,
        maxiter = maxiter,
        tol = tol,
    )

    eof.pcs[:, :] =
        projectfield(eof.eofs[valididx, :], eof.dataset[:, valididx])

    rotatedvar = Statistics.var(eof.pcs, dims = 1)[1, :]
    rotatedvarnorm =
        rotatedvarnorm =
            (rotatedvar ./ sum(rotatedvar)) .* sum(variancefraction(eof))

    sortedidx = sortperm(rotatedvarnorm)[end:-1:1]
    eof.eigenvals[:] = rotatedvar[sortedidx]

    eof.eofs[:, :] = eof.eofs[:, sortedidx]
    eof.pcs[:, :] = eof.pcs[:, sortedidx]

    return
end

"""
Truncate the EOF data structure to trim unneccesay modes
"""
function truncate(eof::EmpiricalOrthogonalFunction, n::Int)
    pcs_trunc = eof.pcs[:, 1:n]
    eofs_trunc = eof.eofs[:, 1:n]
    ev_trunc = eof.eigenvals[1:n]

    return EmpiricalOrthogonalFunction(
        eof.dataset,
        eof.center,
        eof.valididx,
        eofs_trunc,
        ev_trunc,
        pcs_trunc,
    )
end

"""
Apply truncate inplace on EOF data structure
"""
function truncate!(eof::EmpiricalOrthogonalFunction, n::Int)
    eof.pcs = eof.pcs[:, 1:n]
    eof.eofs = eof.eofs[:, 1:n]
    eof.eigenvals = eof.eigenvals[1:n]
    return
end
