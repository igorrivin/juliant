using Primes

struct primePower
    prime::Integer
    power::Int
end
function maxpow(n::Int, p::Int)
    pow = 0
    while true
        rem = n%p
        if rem != 0
            return pow
        end
        pow+=1
        n /=p
    end
end

oddq(n) =  maxpow(n,2)==0

function evenpart(n::Integer)
    ep = 1
    while true
        rem = n%2
        if rem != 0
            return ep
        end
        ep *= 2
        n /= 2
    end
end

function jacobisymbol(a::Integer, b::Integer)
    mp = maxpow(a, 2)
    if mp%2 != 0
        ared =div(a, 2^mp)
        jtemp = jacobisymbol(ared, b)
        b8 = b%8
        if b8 == 1 || b8 == 7
            return jtemp
        else
            return -jtemp
        end
    end
    if a == b-1
        b4 = b%4
        return b4==1?1:3
    end
    if a==1
        return 1
    end
    if a>b
        return jacobisymbol(a%b, b)
    end
    a4 = a%4
    b4 = b%4
    jtemp = jacobisymbol(b, a)
    if a4 == 3 && b4 ==3
         return -jtemp
    else
        return jtemp
    end
end
