# import Base:isapprox, dot

# Below, t_ind returns
# - a tuple if the first argument is a tuple of tuple, and
# - an SVector if the first argument is a tuple of vector.

# Below, earlier methods delegate actions to later methods.  Maybe they need to be
# implemented separately for speed?

# t_ind(t::Tuple3{T}, ind::Tuple3{Int}) where {T} = (t[ind[1]], t[ind[2]], t[ind[3]])
# t_ind(t::Tuple32{T}, i::Int, j::Int, k::Int) where {T} = (t[1][i], t[2][j], t[3][k])  # i, j, k = 1 or 2

# From a tuple of two tuples 1 and 2, each with K entries, construct a tuple (with K entries)
# whose kth entry is the kth entry of either tuple 1 or tuple 2.
@inline t_ind(t::Tuple23, i₁₂::T, j₁₂::T, k₁₂::T) where {T<:Union{GridType,Sign,Integer}} = t_ind(t, (i₁₂,j₁₂,k₁₂))  # i₁₂, j₁₂, k₁₂ = 1 or 2
@inline t_ind(t::Tuple2{NTuple{K}}, ind₁₂::SVector{K,T}) where {K,T<:Union{GridType,Sign,Integer}} = t_ind(t, ind₁₂.data)
@inline t_ind(t::Tuple2{NTuple{K}}, ind₁₂::NTuple{K,T}) where {K,T<:Union{GridType,Sign,Integer}} =  # ind₁₂[k] = 1 or 2
    map((t₁,t₂,i) -> Int(i)==1 ? t₁ : t₂, t[1], t[2], ind₁₂)  # NTuple{K}

# From a tuple of two SVectors 1 and 2, each with K entries, construct an SVector (with K
# entries) whose kth entry is the kth entry of either SVector 1 or SVector 2.
@inline t_ind(t::Tuple2{SVec3}, i₁₂::T, j₁₂::T, k₁₂::T) where {T<:Union{GridType,Sign,Integer}} = t_ind(t, (i₁₂,j₁₂,k₁₂))
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::NTuple{K,T}) where {K,T<:Union{GridType,Sign,Integer}} = t_ind(t, SVector(ind₁₂))
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::SVector{K,T}) where {K,T<:Union{GridType,Sign,Integer}} =  # ind₁₂[k] = 1 or 2
    map((t₁,t₂,i) -> Int(i)==1 ? t₁ : t₂, t[1], t[2], ind₁₂)  # SVector{K}

# From a tuple of K vectors, construct a vector with K entries whose kth entry is the ind[k]-th
# entry of the kth vector.
@inline t_ind(t::Tuple3{AbsVec}, i::Int, j::Int, k::Int) = t_ind(t, (i,j,k))
@inline t_ind(t::Tuple2{AbsVec}, i::Int, j::Int) = t_ind(t, (i,j))
@inline t_ind(t::NTuple{K,AbsVec}, ind::NTuple{K,Int}) where {K} = t_ind(t, SVector(ind))
@inline t_ind(t::NTuple{K,AbsVec}, ind::SVector{K,Int}) where {K} = map((tₖ,iₖ) -> tₖ[iₖ], SVector(t), ind)

# From a tuple of K vectors, construct a tuple of two vectors with K entries whose kth entry
# is taken from the kth vector.
@inline t_ind(t::NTuple{K,AbsVec}, ind::Tuple2{NTuple{K,Int}}) where {K} = (t_ind(t,ind[nN]), t_ind(t,ind[nP]))
@inline t_ind(t::NTuple{K,AbsVec}, ind::Tuple2{SVector{K,Int}}) where {K} = (t_ind(t,ind[nN]), t_ind(t,ind[nP]))

# getindex(t::Tuple3{T}, ind::Tuple3{Int}) where {T} = (t[ind[1]], t[ind[2]], t[ind[3]])
# getindex(t::Tuple32{T}, i::Int, j::Int, k::Int) where {T} = (t[1][i], t[2][j], t[3][k])  # i, j, k = 1 or 2
# getindex(t::Tuple23{T}, i::Int, j::Int, k::Int) where {T} = (t[i][1], t[j][2], t[k][3])  # i, j, k = 1 or 2
# getindex(t::Tuple3{AbsVec{T}}, i::Int, j::Int, k::Int) where {T} = (t[1][i], t[2][j], t[3][k])
# getindex(t::Tuple3{AbsVec{T}}, i::Tuple2{Int}, j::Tuple2{Int}, k::Tuple2{Int}) where {T} =
#     ((t[1][i[1]], t[1][i[2]]), (t[2][j[1]], t[2][j[2]]), (t[3][k[1]], t[3][k[2]]))

# rand3() = (rand(), rand(), rand())
# randn3() = (randn(), randn(), randn())
#
# dot(t1::Tuple3{Real}, t2::Tuple3{Real}) = t1[1]*t2[1] + t1[2]*t2[2] + t1[3]*t2[3]
#
# function norm2(x)
#     if isa(x, Number)
#         return abs(x)
#     else
#         sx = start(x)
#
#         result = 0.0
#         (xi, sx) = next(x, sx)
#         result += norm2(xi)^2
#         while !done(x, sx)
#             (xi, sx) = next(x, sx)
#             result += norm2(xi)^2
#         end
#
#         return sqrt(result)
#     end
# end
#
# function norm2diff(x, y)
#     if isa(x, Number)
#         if !isa(y, Number)
#             return Inf
#         else
#             return abs(x-y)
#         end
#     elseif (isa(x,AbsArr) && isa(y,Tuple)) || (isa(x,Tuple) && isa(y,AbsArr))
#         return Inf
#     else
#         sx, sy = start(x), start(y)
#
#         result = 0.0
#         (xi, sx), (yi, sy) = next(x, sx), next(y, sy)
#         result += norm2diff(xi, yi)^2
#         while !done(x, sx) && !done(y, sy)
#             (xi, sx), (yi, sy) = next(x, sx), next(y, sy)
#             result += norm2diff(xi, yi)^2
#         end
#
#         if !done(x, sx) || !done(y, sy)
#             return Inf
#         else
#             return sqrt(result)
#         end
#     end
# end
#
# isapprox_kernel(x, y, rtol, atol) = norm2diff(x,y) ≤ atol + rtol * max(norm2(x), norm2(y))
# isapprox(x, y; rtol::Real=Base.rtoldefault(Float), atol::Real=eps(Float)) = isapprox_kernel(x, y, rtol, atol)
#
# # Compare x and y with respect to some large number L.  Useful when x or y is zero.
# # It is OK to use ≈ or ≉ for positive numbers (like ∆).
isapprox_wrt(x, y, L::Number, rtol::Real=Base.rtoldefault(Float)) = norm(x-y) ≤ rtol * abs(L)
# isapprox_wrt(x, y, L::Real, rtol::Real=Base.rtoldefault(Float)) = norm2diff(x,y) ≤ rtol * L


function movingavg(l::AbsVec{T}) where {T<:Number}
    # Return (l[1:end-1] + l[2:end]) / 2

    n = length(l)
    if n ≤ 1  # n = 0 or 1
        lmov = float(T)[]
    else  # n ≥ 2
        lmov = Vector{float(T)}(n-1)
        for i = 1:n-1
            lmov[i] = (l[i] + l[i+1]) / 2
        end
    end

    return lmov
end

# # When this is enabled to support isapprox for arrays of tuples, gen_sublprim1d fails
# # at the line where curr[1:2] ≈ prev[end-1:end].  Don't understand why, because
# # the line uses isapprox for arrays of numbers.
# function isapprox(x::AbsArr{NTuple{N,T}}, y::AbsArr{NTuple{N,T}}; rtol::Real=Base.rtoldefault(Float), atol::Real=eps(Float)) where {N,T<:Number}
#     return isapprox_kernel(x, y, rtol, atol)
# end

# function newton(f::Function, f′::Function, x0::Real; rtol::Real=Base.rtoldefault(Float), atol::Real=eps(Float))
#     xₙ = x0
#     fₙ = f(xₙ)
#     while abs(fₙ) > eps(Float) || xₙ ≉ xₙ₋₁  # ensures relative error between xₙ and xₙ₋₁ is less that Base.rtoldefault(Float)
#         xₙ₋₁ = xₙ
#         xₙ -= fₙ / f′(xₙ)
#         fₙ = f(xₙ)
#     end
#
#     return xₙ
# end

# Numerical differentiation by forward difference.
f′fwd(f::Function, xₙ::Number, fₙ::Number=f(xₙ); h::Real=abs(xₙ)*Base.rtoldefault(Float)) = (f(xₙ+h) - fₙ) / h

function newtsol(x₀::Number, f::Function, f′::Function=(x,fₙ)->f′fwd(f,x,fₙ); rtol::Real=Base.rtoldefault(Float), atol::Real=eps(Float))
    # Solve f(x) = 0 using the Newton-Armijo method.
    # x₀: initial guess
    # f: function of x
    # f′: derivative of f at x.  If f′ can be evaluated using the value of f(x),
    #   write f′ such that it takes f(x) as the 2nd argument

    isconverged = true  # true if solution converged; false otherwise
    const maxit = 100    # maximum iterations
    const maxitls = 20   # maximum iterations inside the line search
    const α = 1e-4
    const maxs = 1 / Base.rtoldefault(Float)  # maximum step size
    const perturbls = eps(Float)^0.75  # ≈ 1e-12, between sqrt(eps) and eps; some perturbation allowed in line search

    # Initialize.
    n = 0
    xₙ = float(x₀)
    const T = typeof(xₙ)
    fₙ::T = f(xₙ)
    const τf = rtol*abs(fₙ) + atol
    const rx₀ = abs(x₀)
    const has2ndarg = method_exists(f′, Tuple2{Number})

    # Perform the Newton method.
    # info("fₙ = $fₙ")
    while abs(fₙ) > τf
        # info("n = $n")

        # abs(xₙ/x₀) ≤ 1e3 || (isconverged = false; break)
        # abs(xₙ/x₀) ≤ 1e3 || throw(ErrorException("Solution xₙ = $xₙ has diverged from x₀ = $x₀."))

        λ = 1.
        nls = 0  # line search iteration counter

        f′ₙ::T = has2ndarg ? f′(xₙ,fₙ) : f′(xₙ)

        # Avoid too large Newton steps.
        s = -fₙ/f′ₙ
        # info("fₙ = $fₙ, f′ₙ = $f′ₙ, xₙ = $xₙ, s = $s")
        abs(s) ≤ maxs || (isconverged = false; break)
        # abs(s) ≤ maxs || throw(ErrorException("Newton step s = $s is larger than maximum step size $maxs."))
        xₙ₊₁ = xₙ + λ*s
        fₙ₊₁ = f(xₙ₊₁)

        # Perform the line search to determine λ.  The stopping criterion does not
        # have perturbls ≈ 1e-12 on the RHS, but I guess this kind of perturbation
        # allows update in xₙ even in the situation where line search is supposed
        # to fail.
        while abs(fₙ₊₁) ≥ (1 - α*λ) * abs(fₙ) + perturbls
            # info("nls = $nls, fₙ₊₁ = $(fₙ₊₁)")
            λ /= 2
            xₙ₊₁ = xₙ + λ*s
            fₙ₊₁::T = f(xₙ₊₁)
            nls += 1

            # Too many iteration steps in line search
            nls ≤ maxitls || (isconverged = false; break)
            # nls ≤ maxitls || throw(ErrorException("Line search fails in $nls iteration steps."))
        end

        # Step accepted; continue the Newton method.
        xₙ = xₙ₊₁
        n += 1

        # Too many iteration steps in Newton's method.
        n ≤ maxit || (isconverged = false; break)
        # n ≤ maxit || throw(ErrorException("Newton method fails to converge in $n iteration steps."))

        fₙ = f(xₙ)
        # info("fₙ = $fₙ")
    end
    # info("Newton done")

    return xₙ, isconverged
end
