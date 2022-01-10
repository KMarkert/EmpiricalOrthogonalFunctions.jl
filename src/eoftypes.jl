# This file defines the EmpiricalOrthogonalFunction types

# abstract EOF type
"""
Abstract type to build subtypes from
"""
abstract type AbstractEmpiricalOrthogonalFunction end

# eof struct for 2d data (space, time)
"""
    EmpiricalOrthogonalFunction(dataset; center=true, ddof=1)

Create an Empirical Orthogonal Function object.
The EOF solution is computed at initialization time. Method
calls are used to retrieve or update computed quantities.
"""
mutable struct EmpiricalOrthogonalFunction <: AbstractEmpiricalOrthogonalFunction
    ""
    dataset::AbstractArray{<:Any,2}

    ""
    center::AbstractArray{<:Any,2}

    valididx::AbstractArray{Int64,1}

    eofs::AbstractArray{<:Any,2}

    eigenvals::AbstractArray{<:Any,1}

    pcs::AbstractArray{<:Any,2}

    EmpiricalOrthogonalFunction(
        dataset::AbstractArray{<:Any,2},
        center::AbstractArray{<:Any,2},
        valididx::AbstractArray{Int64,1},
        eofs::AbstractArray{<:Any,2},
        eigenvals::AbstractArray{<:Any,1},
        pcs::AbstractArray{<:Any,2},
    ) = new(dataset, center, valididx, eofs, eigenvals, pcs)

    """
    """
    function EmpiricalOrthogonalFunction(
        dataset::AbstractArray{<:Any,2};
        center::Bool = true,
        ddof::Int = 1,
    )

        # Store information about the shape/size of the input data.
        shape = size(dataset)
        timelen = shape[1]

        if shape[1] > shape[2]
            @error(
                "time dimension has more records than the space dimensions and an EOF cannot be calculated, time = $(shape[1]), space= $(shape[2])"
            )
            return
        end

        # Remove the time mean of the input data unless explicitly told
        # not to by the "center" argument.
        if center
            spatial_mean = Statistics.mean(dataset, dims = 1)
            centered = @. dataset - spatial_mean
        else
            spatial_mean = zeros(shape[1])
            centered = dataset
        end

        # Find the indices of values that are not missing in one row. All the
        # rows will have missing values in the same places provided the
        # array was centered. If it wasn't then it is possible that some
        # missing values will be missed and the singular value decomposition
        # will produce not a number for everything.
        isnotmissing = findall(x -> x !== missing, spatial_mean[1, :])

        # Compute the singular value decomposition of the design matrix.
        U, S, V = svd(centered[:, isnotmissing])

        # Singular values are the square-root of the eigenvalues of the
        # covariance matrix. Construct the eigenvalues appropriately and
        # normalize by N-ddof where N is the number of observations. This
        # corresponds to the eigenvalues of the normalized covariance matrix.
        norm = timelen - 1
        eigenvals = @. S * S / norm

        eofs = Array{Union{Missing,Float64}}(missing, shape[2], shape[1])
        eofs[isnotmissing, :] = V[:, :]

        # Remove the scaling on the principal component time-series that is
        # implicitily introduced by using SVD instead of eigen-decomposition.
        # The PCs may be re-scaled later if required.
        pcs = similar(U)
        for i = 1:length(S)
            pcs[:, i] = U[:, i] * S[i]
        end

        EmpiricalOrthogonalFunction(
            dataset,
            spatial_mean,
            isnotmissing,
            eofs,
            eigenvals,
            pcs,
        )
    end

    """

    """
    function EmpiricalOrthogonalFunction(
        dataset::AbstractArray{<:Any,3};
        timedim::Int = 3,
        kwargs...,
    )
        geodims = [x for x in 1:3 if x != timedim]
        shape = size(dataset)
        spatial_shape = shape[geodims]
        time_shape = shape[timedim]

        if timedim != 1
            dataset = permutedims(dataset, (timedim, geodims...))
        end

        flattened = reshape(dataset, (time_shape, prod(spatial_shape)))

        EmpiricalOrthogonalFunction(flattened; kwargs...)
    end

end
