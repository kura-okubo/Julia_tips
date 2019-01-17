module Tools
export snr, smooth, nextpow2, abs_max, standardize, mad, savitsky_golay,
     running_mean, pws

using ..SeisJul

"""
    snr(A,sampling_rate)

Signal to noise ratio of cross-correlations in matrix `A` with `sampling_rate'.

Follows method of Clarke et. al, 2011. Measures SNR at each point.
"""
function snr(A::AbstractArray,sampling_rate::Real)
    t,N = size(A)
    A_mean = mean(A,dims=2)

    # calculate noise and envelope functions
    sigma = mean(A.^2,dims=2) .- A_mean.^2
    sigma = sqrt.(sigma./ (N-1))
    return sigma
end

"""
    smooth(x,half_win=3)

Smooth vector `x` with half-window `half_win` (defaults to 3).
"""
function smooth(x::AbstractArray; half_win::Int=3)
    # only use odd numbers
    if half_win % 2 != 1
        half_win += oneunit(half_win)
    end
    window_len = 2 * half_win + 1

    # extending the data at beginning and at the end
    # to apply the window at the borders
    A = vcat(x,x[window_len:-1:2])
    append!(A,x[end-1:-1:end-window_len+1])

    # convolve with boxcar
    w = rect(window_len)
    y = conv(A,w/sum(w))
    y = y[window_len+half_win:end-window_len-half_win+1]
end

"""
    nextpow2(x)

Returns the next power of 2 of real, positive number `x`.
"""
function nextpow2(x::Real)
    ceil(Int,log2(abs(x)))
end

"""
    abs_max(A)

Returns array `A` divided by its absolute maximum value.
"""
function abs_max(A::AbstractArray)
    A ./ maximum(abs.(A),dims=1)
end

"""
    standardize(A)

Demean and standardize array `A` to unit std.
"""
function standardize(A::AbstractArray)
    (A .- mean(A,dims=1)) ./ std(A,dims=1)
end

"""
    mad(A)

Median Absolute Deviation of array `A`.
MAD = median(|Xi- median(A)|)
"""
function mad(A::AbstractArray)
    median(abs.(A .- median(A,dims=1)),dims=1)

end

"""
    savitsky_golay(x, window, polyOrder, [deriv])

Polynomial smoothing of vector `x` with Savitsky Golay filter.
Polynomial order `polyOrder' must be less than window length `N`.
Modified from https://github.com/BBN-Q/Qlab.jl/blob/master/src/SavitskyGolay.jl
"""
function savitsky_golay(x::Vector, N::Int, polyOrder::Int; deriv::Int=0)
    #Some error checking
    @assert isodd(N) "Window size must be an odd integer."
    @assert polyOrder < N "Polynomial order must me less than window size."

    halfWindow = Int((N-1)/2)

    #Setup the S matrix of basis vectors.
    S = zeros(N, polyOrder+1)
    for ct = 0:polyOrder
        S[:,ct+1] = Array(-halfWindow:halfWindow).^ct
    end

    #Compute the filter coefficients for all orders
    G = S*pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv);

    #Pad the signal with the endpoints and convolve with filter
    paddedX = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
    y = conv(filterCoeffs[end:-1:1], paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]
end

"""
    running_mean(x,N)

Returns array `x` smoothed by running mean of `N` points.
If N is even, reduces N by 1.
"""
function running_mean(x::AbstractArray,N::Int)
    if N % 2 == 0
        N -= 1
    end
    halfWindow = div(N,2)
    paddedx = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
    y = conv(paddedx,ones(N) / N)
    return  y[2*halfWindow+1:end-2*halfWindow]
end

"""
    pws(A,sampling_rate,[power],[timegate])

Performs phase-weighted stack on array `A` of time series.

Follows methods of Schimmel and Paulssen, 1997.
If s(t) is time series data,
S(t) = s(t) + i*H(s(t)), where H(s(t)) is Hilbert transform of s(t)
S(t) = s(t) + i*H(s(t)) = A(t)*exp(i*phi(t)), where
A(t) is envelope of s(t) and phi(t) is phase of s(t)
Phase-weighted stack, g(t), is then:
g(t) = 1/N sum j = 1:N s_j(t) * | 1/N sum k = 1:N exp[i * phi_k(t)]|^v
where N is number of traces used, v is sharpness of phase-weighted stack

"""
function pws(A::AbstractArray, sampling_rate::Real=20; power::Int=2, timegate::Int=5)
    M,N = size(A)
    analytic = A .+ im .* hilbert(A)
    phase = angle.(analytic)
    phase_stack = mean(exp.(im.*phase),dims=2)[:,1] ./ N # reduce array dimension
    phase_stack = abs.(phase_stack).^2

    # smoothing
    timegate_samples = Int(timegate * sampling_rate)
    phase_stack = running_mean(phase_stack,timegate_samples)
    weighted = A * phase_stack
    return weighted
end

end
