# # WAVELET  1D Wavelet transform with optional significance testing
#   wave, period, scale, coi = wavelet(Y, dt, pad, dj, s0, J1, mother, param)
#
#   Computes the wavelet transform of the vector Y (length N),
#   with sampling rate DT.
#
#   By default, the Morlet wavelet (k0=6) is used.
#   The wavelet basis is normalized to have total energy=1 at all scales.
#
# INPUTS:
#
#    Y = the time series of length N.
#    DT = amount of time between each Y value, i.e. the sampling time.
#
# OUTPUTS:
#
#    WAVE is the WAVELET transform of Y. This is a complex array
#    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
#    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
#    The WAVELET power spectrum is ABS(WAVE)**2.
#    Its units are sigma**2 (the time series variance).
#
# OPTIONAL INPUTS:
#
# *** Note *** if none of the optional variables is set up, then the program
#   uses default values of -1.
#
#    PAD = if set to 1 (default is 0), pad time series with zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#
#    DJ = the spacing between discrete scales. Default is 0.25.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    S0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    J1 = the # of scales minus one. Scales range from S0 up to S0*2**(J1*DJ),
#        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
#
#    MOTHER = the mother wavelet function.
#             The choices are 'MORLET', 'PAUL', or 'DOG'
#
#    PARAM = the mother wavelet parameter.
#            For 'MORLET' this is k0 (wavenumber), default is 6.
#            For 'PAUL' this is m (order), default is 4.
#            For 'DOG' this is m (m-th derivative), default is 2.
#
#
# OPTIONAL OUTPUTS:
#
#    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
#           to the SCALEs.
#
#    SCALE = the vector of scale indices, given by S0*2**(j*DJ), j=0...J1
#            where J1+1 is the total # of scales.
#
#    COI = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.

using FFTW
using SpecialFunctions
using Optim
using Statistics

function wavelet(Y::Vector{Float64}, dt::Float64; pad=0, dj=-1, s0=-1, J1=-1, mother=-1, param=-1, freq=nothing)
    n1 = length(Y)
    if s0 == -1
        s0 = 2 * dt
    end
    if dj == -1
        dj = 1. / 4.
    end
    if J1 == -1
        J1 = trunc((log(n1 * dt / s0) / log(2)) / dj)
    end
    if mother == -1
        mother = "MORLET"
    end
    
    # construct time series to analyze, pad if necessary
    x = Y .- mean(Y)
    if pad == 1 
        # power of 2 nearest to N
        base2 = round(log(n1) / log(2) + 0.4999)
        nzeroes = Int64(2 ^ (base2 + 1) - n1)
        x = vcat(x, zeros(nzeroes))
    end
    n = length(x)
    
    # construct wavenumber array used in transform [Eqn(5)]
    kplus = collect(1:trunc(n / 2))
    kplus = kplus .* (2pi / (n * dt))
    kminus = collect(1:trunc((n - 1) / 2))
    kminus = sort(-kminus .* (2pi / (n * dt)))
    k = vcat([0.0], kplus, kminus)
    
    # compute FFT of the (padded) time series
    f = fft(x)  # [Eqn(3)]
    
    # construct SCALE array & empty PERIOD & WAVE arrays
    if uppercase(mother) == "MORLET"
        if param == -1
            param = 6
        end
        fourier_factor = 4 * pi / (param + sqrt(2 + param^2))
    elseif uppercase(mother) == "PAUL"
        if param == -1
            param = 4
        end
        fourier_factor = 4 * pi / (2 * param + 1)
    elseif uppercase(mother) == "DOG"
        if param == -1
            param = 2
        end
        fourier_factor = 2pi * sqrt(2 / (2 * param + 1))
    else
        fourier_factor = NaN
    end
    if isnothing(freq)
        j = collect(0:J1)
        scale = s0 .* 2.0 .^(j .* dj)
        freq = 1 ./ (fourier_factor .* scale)
        period = 1 ./ freq
    else
        scale = 1 ./ (fourier_factor .* freq)
        period = 1 ./ freq
    end
    # define the wavelet array
    
    wave = zeros(Complex{Float64}, length(scale), n)
    
    # loop through all scales and compute transform
    for a1 in eachindex(scale)
        daughter, fourier_factor, coi, _ = wave_bases(mother, k, scale[a1], param)
        wave[a1, :] = ifft(f .* daughter)  # wavelet transform[Eqn(4)]
    end
    
    # COI [Sec.3g]
    coi = coi * dt .* vcat(insert!(collect(0:trunc((n1 + 1) / 2) - 2), 1, 1E-5), insert!(reverse(collect(0:trunc(n1 / 2) - 2)),length(reverse(collect(0:trunc(n1 / 2) - 2))), 1E-5)) #push!(reverse(collect(0:trunc(n1 / 2) - 2)),1E-5)
    wave = wave[:, 1:n1]  # get rid of padding before returning
    return wave, period, scale, coi
end

# --------------------------------------------------------------------------
# WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
#
#  DAUGHTER,FOURIER_FACTOR,COI,DOFMIN = wave_bases(MOTHER,K,SCALE,PARAM)
#
#   Computes the wavelet function as a function of Fourier frequency,
#   used for the wavelet transform in Fourier space.
#   (This program is called automatically by WAVELET)
#
# INPUTS:
#
#    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
#    K = a vector, the Fourier frequencies at which to calculate the wavelet
#    SCALE = a number, the wavelet scale
#    PARAM = the nondimensional parameter for the wavelet function
#
# OUTPUTS:
#
#    DAUGHTER = a vector, the wavelet function
#    FOURIER_FACTOR = the ratio of Fourier period to scale
#    COI = a number, the cone-of-influence size at the scale
#    DOFMIN = a number, degrees of freedom for each point in the wavelet power
#             (either 2 for Morlet and Paul, or 1 for the DOG)

function wave_bases(mother, k, scale, param)
    n = length(k)
    kplus = k .> 0.0
    if mother == "MORLET"  # -----------------------------------  Morlet
        if param == -1
            param = 6
        end
        k0 = copy(param)
        # calc psi_0(s omega) from Table 1
        expnt = -(scale .* k .- k0) .^ 2 / 2.0 .* kplus
        norm = sqrt(scale * k[2]) * (pi ^ (-0.25)) * sqrt(n)
        daughter = norm * exp.(expnt)
        daughter = daughter .* kplus  # Heaviside step function
        # Scale-->Fourier [Sec.3h]
        fourier_factor = 4pi / (k0 + sqrt(2 + k0 ^ 2))
        coi = fourier_factor / sqrt(2)  # Cone-of-influence [Sec.3g]
        dofmin = 2  # Degrees of freedom
    elseif mother == "PAUL"  # --------------------------------  Paul
        if param == -1
            param = 4
        end
        m = param
        # calc psi_0(s omega) from Table 1
        expnt = -scale .* k .* kplus
        norm_bottom = sqrt(m * prod(collect(1:(2 * m)-1)))
        norm = sqrt(scale * k[2]) * (2 ^ m / norm_bottom) * sqrt(n)
        daughter = norm .* ((scale .* k) .^ m) .* exp.(expnt) .* kplus
        fourier_factor = 4pi / (2 * m + 1)
        coi = fourier_factor * sqrt(2)
        dofmin = 2
    elseif mother == "DOG"  # --------------------------------  DOG
        if param == -1
            param = 2
        end
        m = param
        # calc psi_0(s omega) from Table 1
        expnt = -(scale * k) .^ 2 ./ 2.0
        norm = sqrt(scale * k[2] / gamma(m + 0.5)) * sqrt(n)
        daughter = -norm * (1im ^ m) * ((scale .* k) .^ m) .* exp.(expnt)
        fourier_factor = 2pi * sqrt(2. / (2 * m + 1))
        coi = fourier_factor / sqrt(2)
        dofmin = 1
    else
        println("Mother must be one of MORLET, PAUL, DOG")
    end
    return daughter, fourier_factor, coi, dofmin
end

# --------------------------------------------------------------------------
# WAVE_SIGNIF  Significance testing for the 1D Wavelet transform WAVELET
#
#   SIGNIF = wave_signif(Y,DT,SCALE,SIGTEST,LAG1,SIGLVL,DOF,MOTHER,PARAM)
#
# INPUTS:
#
#    Y = the time series, or, the VARIANCE of the time series.
#        (If this is a single number, it is assumed to be the variance...)
#    DT = amount of time between each Y value, i.e. the sampling time.
#    SCALE = the vector of scale indices, from previous call to WAVELET.
#
#
# OUTPUTS:
#
#    SIGNIF = significance levels as a function of SCALE
#    FFT_THEOR = output theoretical red-noise spectrum as fn of PERIOD
#
#
# OPTIONAL INPUTS:
#    SIGTEST = 0, 1, or 2.    If omitted, then assume 0.
#
#         If 0 (the default), then just do a regular chi-square test,
#             i.e. Eqn (18) from Torrence & Compo.
#         If 1, then do a "time-average" test, i.e. Eqn (23).
#             In this case, DOF should be set to NA, the number
#             of local wavelet spectra that were averaged together.
#             For the Global Wavelet Spectrum, this would be NA=N,
#             where N is the number of points in your time series.
#         If 2, then do a "scale-average" test, i.e. Eqns (25)-(28).
#             In this case, DOF should be set to a
#             two-element vector [S1,S2], which gives the scale
#             range that was averaged together.
#             e.g. if one scale-averaged scales between 2 and 8,
#             then DOF=[2,8].
#
#    LAG1 = LAG 1 Autocorrelation, used for SIGNIF levels. Default is 0.0
#
#    SIGLVL = significance level to use. Default is 0.95
#
#    DOF = degrees-of-freedom for signif test.
#         IF SIGTEST=0, then (automatically) DOF = 2 (or 1 for MOTHER='DOG')
#         IF SIGTEST=1, then DOF = NA, the number of times averaged together.
#         IF SIGTEST=2, then DOF = [S1,S2], the range of scales averaged.
#
#       Note: IF SIGTEST=1, then DOF can be a vector (same length as SCALEs),
#            in which case NA is assumed to vary with SCALE.
#            This allows one to average different numbers of times
#            together at different scales, or to take into account
#            things like the Cone of Influence.
#            See discussion following Eqn (23) in Torrence & Compo.
#
#    GWS = global wavelet spectrum, a vector of the same length as scale.
#          If input then this is used as the theoretical background spectrum,
#          rather than white or red noise.

function wave_signif(variance, dt, scale::Vector{Float64}; sigtest=0, lag1=0.0, siglvl=0.95,
                dof=[nothing], mother="MORLET", param=nothing, gws=nothing)
#function wave_signif(variance, dt, scale::Vector{Float64}; sigtest=0, lag1=0.0, siglvl=0.95,
#                dof=[nothing], mother="MORLET", param=nothing, gws=nothing)
    J1 = length(scale) - 1
    dj = log2(scale[2] / scale[1])

    # get the appropriate parameters [see Table(2)]
    if mother == "MORLET"  # ----------------------------------  Morlet
        empir = [2.0, -1.0, -1.0, -1.0]
        if isnothing(param)
            param = 6
            empir[2:end] = [0.776, 2.32, 0.60]
        end
        k0 = param
        # Scale-->Fourier [Sec.3h]
        fourier_factor = 4pi / (k0 + sqrt(2 + k0 ^ 2))
    elseif mother == "PAUL"
        empir = [2.0, -1.0, -1.0, -1.0]
        if isnothing(param)
            param = 4
            empir[2:end] = [1.132, 1.17, 1.5]
        end
        m = param
        fourier_factor = 4pi / (2 * m + 1)
    elseif mother == "DOG"  # -------------------------------------Paul
        empir = [1.0, -1.0, -1.0, -1.0]
        if isnothing(param)
            param = 2
            empir[2:end] = [3.541, 1.43, 1.4]
        elseif param == 6  # --------------------------------------DOG
            empir[2:end] = [1.966, 1.37, 0.97]
        end
        m = param
        fourier_factor = 2pi * sqrt(2. / (2 * m + 1))
    else
        println("Mother must be one of MORLET, PAUL, DOG")
    end
    
    period = scale * fourier_factor
    dofmin = empir[1]  # Degrees of freedom with no smoothing
    Cdelta = empir[2]  # reconstruction factor
    gamma_fac = empir[3]  # time-decorrelation factor
    dj0 = empir[4]  # scale-decorrelation factor
    freq = dt ./ period  # normalized frequency
    if !isnothing(gws) # use global-wavelet as background spectrum
        fft_theor = gws
    else
        # [Eqn(16)]
        fft_theor = (1 - lag1 ^ 2) ./ (1 .- 2lag1 .* cos.(2pi * freq) .+ lag1 ^ 2)
        fft_theor = variance * fft_theor  # include time-series variance

    end
    signif = fft_theor
    if isnothing(dof[begin])
        dof=dofmin   #dof = [dofmin]
    end
    if sigtest == 0  # no smoothing, DOF=dofmin [Sec.4]
        dof = dofmin
        chisquare = chisquare_inv(siglvl, dof) / dof
        signif = fft_theor * chisquare  # [Eqn(18)]
    elseif sigtest == 1  # time-averaged significance
        if length(dof) == 1
            dof = zeros(J1) .+ dof
        end
        dof[dof .< 1] .= 1
        # [Eqn(23)]
        dof = dofmin .* sqrt.(1 .+ (dof .* dt ./ gamma_fac ./ scale) .^ 2)
        dof[dof .< dofmin] .= dofmin   # minimum DOF is dofmin
        for a1 in 1:J1 + 1
            chisquare = chisquare_inv(siglvl, dof[a1]) / dof[a1]
            signif[a1] = fft_theor[a1] * chisquare
        end
    elseif sigtest == 2  # time-averaged significance
        if length(dof) != 2
            println("ERROR: DOF must be set to [S1,S2], the range of scale-averages")
        end
        if Cdelta == -1
            println("ERROR: Cdelta & dj0 not defined for $mother with param = $param")
        end
        s1 = dof[1]
        s2 = dof[2]
        avg = (scale .>= 2) .& (scale .< 8)  # scales between S1 & S2
        navg = sum(collect((scale .>= 2) .& (scale .< 8)))
        if navg == 0
            println("ERROR: No valid scales between $s1 and $s2")
        end
        Savg = 1. / sum(1. ./ scale[avg])  # [Eqn(25)]
        Smid = exp((log(s1) + log(s2)) / 2.)  # power-of-two midpoint
        dof = (dofmin * navg * Savg / Smid) * sqrt(1 + (navg * dj / dj0) ^ 2)  # [Eqn(28)]
        fft_theor = Savg * sum(fft_theor[avg] ./ scale[avg])  # [Eqn(27)]
        chisquare = chisquare_inv(siglvl, dof) / dof
        signif = (dj * dt / Cdelta / Savg) * fft_theor * chisquare   # [Eqn(26)]
    else
        println("ERROR: sigtest must be either 0, 1, or 2")
    end
    return signif
end

# --------------------------------------------------------------------------
# CHISQUARE_INV  Inverse of chi-square cumulative distribution function (cdf).
#
#   X = chisquare_inv(P,V) returns the inverse of chi-square cdf with V
#   degrees of freedom at fraction P.
#   This means that P*100 percent of the distribution lies between 0 and X.
#
#   To check, the answer should satisfy:   P==gammainc(X/2,V/2)

# Uses FMIN and CHISQUARE_SOLVE

function chisquare_inv(P, V)
    if (1 - P) < 1E-4
        println("P must be < 0.9999")
    end
    if P == 0.95 && V == 2  # this is a no-brainer
        X = 5.9915
        return X
    end
    MINN = 0.01  # hopefully this is small enough
    MAXX = 1  # actually starts at 10 (see while loop below)
    X = 1
    TOLERANCE = 1E-4  # this should be accurate enough

    while (X + TOLERANCE) >= MAXX  # should only need to loop thru once
        MAXX = MAXX * 10.
    # this calculates value for X, NORMALIZED by V    
        chisquare_solveT(XGUESS) = chisquare_solve(XGUESS, P, V)
        res = optimize(chisquare_solveT, MINN, MAXX, rel_tol=TOLERANCE)
        X = Optim.minimizer(res)  #fminbound(chisquare_solve, MINN, MAXX, args=(P, V), xtol=TOLERANCE)
        MINN = MAXX
    end
    X = X * V  # put back in the goofy V factor
    return X  # end of code
end

# --------------------------------------------------------------------------
# CHISQUARE_SOLVE  Internal function used by CHISQUARE_INV
    #
    #   PDIFF=chisquare_solve(XGUESS,P,V)  Given XGUESS, a percentile P,
    #   and degrees-of-freedom V, return the difference between
    #   calculated percentile and P.

    # Uses GAMMAINC
    #
    # Written January 1998 by C. Torrence

    # extra factor of V is necessary because X is Normalized

function chisquare_solve(XGUESS, P, V)
    PGUESS, temp = gamma_inc(V / 2, V * XGUESS / 2)  # incomplete Gamma function #gammainc(V / 2, V * XGUESS / 2)  
    PDIFF = abs(PGUESS - P)  # error in calculated P          
    TOL = 1E-4
    
    if PGUESS >= 1 - TOL  # if P is very close to 1 (i.e. a bad guess)
        PDIFF = XGUESS  # then just assign some big number like XGUESS

    end
    return PDIFF
end


