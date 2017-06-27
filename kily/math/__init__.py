import cmath as math
import random

def prime_test_sieve(a: int) -> bool:
    return all(a % b != 0 for b in range(2, int(math.sqrt(a)+1), 2))

def prime_test_miller_rabin(a: int, iters: int):
    if a == 2:
        return True
    if not a & 1:
        return False

    def check(a, s, d, n):
        x = pow(a, d, n)
        if x == 1:
            return True
        for i in range(s - 1):
            if x == n - 1:
                return True
            x = pow(x, 2, n)
        return x == n - 1

    s = 0
    d = a - 1

    while d % 2 == 0:
        d >>= 1
        s += 1

    for i in range(iters):
        a = random.randrange(2, a - 1)
        if not check(a, s, d, a):
            return False
    return True

def prime_test(a: int, slow_test_cutoff: int = 4*10**6, pow_probability: int = 16):
    if a > slow_test_cutoff:
        real_iters = pow_probability//2
        return prime_test_miller_rabin(a, real_iters)
    else:
        return prime_test_sieve(a)


# return (g, x, y) a*x + b*y = gcd(x, y)
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, x, y = egcd(b % a, a)
        return (g, y - (b // a) * x, x)

def gcd(a: int, b: int) -> int:
    return egcd(min(a, b), max(a, b))[0]


def modulo_inverse(a: int, mod: int) -> int:
    if a > mod:
        a %= mod
    g, inv = egcd(a, mod)
    if g == 1:
        return inv
    else:
        raise ValueError(f"Values of a ({a}) and mod ({mod}) are not co-prime.")