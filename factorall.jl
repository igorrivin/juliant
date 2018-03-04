using Primes
using OffsetArrays

struct primePower
    prime::Int
    power::Int
end

function maxpow(n::Int, p::Int)
    pow = 1
    n /= p
    while true
        rem = n%p
        if rem != 0
            return pow
        end
        pow+=1
        n /=p
    end
end

function factorall(n)
    dictlist = [Dict{Int, Int}() for i in 1:n]
    theprimes = primes(Int(ceil(sqrt(n))))
    for i in theprimes
        for j in i:i:n
            @inbounds dictlist[j][i]=maxpow(j, i)
        end
    end
    for i in 2:n
        dl = dictlist[i]
        if isempty(dl)
            dl[i] = 1
        end
        pfacs = prod(i^dl[i] for i in keys(dl))
        therat = i/pfacs
        if therat > 1
            dl[therat]=1
        end
    end
    return dictlist
end

function factorall(m,n)
    dictlist = [Dict{Int, Int}() for i in m:n]
    dictlist = OffsetArray(dictlist, m:n)
    theprimes = primes(Int(ceil(sqrt(n))))
    for i in theprimes
        begprime = Int(ceil(m/i))*i
        for j in begprime:i:n
            @inbounds dictlist[j][i]=maxpow(j, i)
        end
    end
    thetwo = max(2, m)
    for i in thetwo:n
        dl = dictlist[i]
        if isempty(dl)
            dl[i] = 1
        end
        pfacs = prod(i^dl[i] for i in keys(dl))
        therat = i/pfacs
        if therat > 1
            dl[therat]=1
        end
    end
    return dictlist
end

newfact = memoizeany(factorall)

function domult(ndict, func)
    if isempty(ndict)
        return 1
    end
    return prod(func(i, ndict[i]) for i in keys(ndict))
end

function doallmult(n, func)
    return [domult(i, func) for i in newfact(n)]
end

function doallmult(m, n, func)
    tmpar = collect(factorall(m, n)) #compiler bug
    res = [domult(i, func) for i in tmpar]
    return OffsetArray(res, m:n)
end


function makesigmab(k)
    thedict = Dict((0, 0)=>0)
    return (p, n) -> sigkb(p, n, k, thedict)
end

function makesigma(k)
    thedict = Dict((0, 0)=>0)
    if k==0
        return (p, n)-> n+1
    end
    return (p, n) -> sigk(p, n, k, thedict)
end

function sigkb(p, n, k,thedict)
    if haskey(thedict, (p, n))
        return thedict[(p, n)]
    end
    bp = BigInt(p)
    @inbounds tmp = sum(bp^(i*k) for i in 0:n)
    thedict[(p, n)] = tmp
    return tmp
end

function sigk(p, n, k,thedict)
    if haskey(thedict, (p, n))
        return thedict[(p, n)]
    end
    @inbounds tmp = sum(p^(i*k) for i in 0:n)
    thedict[(p, n)] = tmp
    return tmp
end

function memoizeany(func)
    thedict = Dict()
    return (a)-> memoizeanyaux(a, func, thedict)
end

function memoizeanyaux(a, func, thedict)
    if haskey(thedict, a)
        return thedict[a]
    end
    res = func(a)
    thedict[a] =res
    return res
end

function memoize( func)
    thedict = Dict((0, 0)=>0)
    return (p, n)-> memoizeaux(p, n, func, thedict)
end

function memoizeaux(p, n, func, thedict)
    if haskey(thedict, (p, n))
        return thedict[(p, n)]
    end
    res = func(p, n)
    thedict[(p, n)] = res
    return res
end

mobius(p, n) = n==1?-1:0

phifunc(p, n) = p^(n-1)*(p-1)
