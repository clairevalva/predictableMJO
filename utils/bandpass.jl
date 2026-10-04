# FIR Lanczos filtering along dimension 1 (time). Callers with longitude × time
# arrays must transpose before and after filtering. Final padding and delay
# compensation preserve the original filter implementation.
using DSP

function makeFilterMJO(;days_low = 100, days_high = 20, sample_fs = 1, lanN = 49)
    # make bandpass filter for MJO as in Kessler 2001
    # call as F = makeFilterMJO(), filt(F, x), where the first index of x corresponds to time

    lowband = 1 / days_low 
    highband = 1 / days_high

    responsetype = Bandpass(lowband, highband)

    designmethod = FIRWindow(lanczos(lanN))

    
    F = digitalfilter(responsetype, designmethod; fs = sample_fs)
    return F
end

function makeMJOFF(;days_low = 96, days_high = 20, sample_fs = 1, lanN = 49)

    nPad = round(Int, lanN / 2)
    nTaps = 2*nPad

    F = makeFilterMJO(;days_low = days_low, days_high = days_high, sample_fs = sample_fs, lanN = nTaps)

    function applyF(X)
        # X has first index as time, i.e. size = nT × nD1 × nD2 …
        sizeX = size(X)
        nT = sizeX[1]

        sizeX_pad = zeros(Int, length(sizeX))
        sizeX_pad .= sizeX
        sizeX_pad[1] = nPad

        Xadd = zeros(sizeX_pad...)

        X_pad = vcat(X, Xadd)

        X_filt = filt(F, X_pad)
        Xreturn = selectdim(X_filt , 1, (nPad + 1):(nT + nPad))

        return Xreturn
    end

    return applyF
end
