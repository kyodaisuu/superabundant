#!/usr/bin/env python3
"""
Superabundant number search using the structural theorems of
Alaoglu and Erdos (1944), "On Highly Composite and Similar Numbers",
Trans. Amer. Math. Soc. 56, 448-469.

A number n is superabundant if sigma(m)/m < sigma(n)/n for all 0 < m < n,
where sigma(n) is the sum of divisors of n.

Theorem 1 (Alaoglu-Erdos): if n = 2^k2 * 3^k3 * 5^k5 * ... * p^kp is
superabundant, then k2 >= k3 >= k5 >= ... >= kp. In particular a
superabundant number only ever uses the first k primes (2, 3, 5, ..., p_k)
with no gaps, and the exponents are non-increasing as the prime grows.

Theorem 2 (Alaoglu-Erdos): let q < r be primes, both possibly dividing a
superabundant number n, with k_q the exponent of q in n. Set
beta = floor(k_q * log(q) / log(r)). Then the exponent k_r of r in n is
one of beta - 1, beta, beta + 1.

Theorem 3 (Alaoglu-Erdos): if p is the largest prime factor of a
superabundant number n, then the exponent of p is 1, except for n = 4
and n = 36.

Instead of testing every integer up to a bound N (which is only feasible
for N as large as roughly 10^16 with a naive sieve-and-factor approach),
we use Theorem 1 to restrict the search to non-increasing exponent
sequences (a_1 >= a_2 >= ... >= a_k >= 1) applied to the first k primes,
Theorem 2 to narrow, at every step, the exponent of each new prime down
to at most three candidate values (instead of the full range allowed by
Theorem 1 alone), and Theorem 3 to discard, along the way, candidates
whose top exponent already violates the theorem. Theorem 2 is applied
with q fixed as the prime 2: once a1 (the exponent of 2) is chosen at
the top of the search, it anchors a tight predicted exponent for every
later prime, which turns what would otherwise be a search branching
factor of "up to the previous exponent" into a branching factor of "at
most 3" at every step. This is the same reduction that let Alaoglu and
Erdos build their own table of superabundant numbers up to 10^18 by
hand in 1944, and it lets this program comfortably reach far larger
bounds than Theorem 1 alone would allow, while still working with
arbitrary-precision Python integers so N can be as large as 10^n for
any n.

Usage:
    python superabundant.py 16      # search up to 10^16 (default: 16)
"""

import sys


def sieve_primes(limit):
    """Return the list of primes <= limit using a simple sieve of Eratosthenes."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = bytearray(len(range(i * i, limit + 1, i)))
    return [i for i, flag in enumerate(is_prime) if flag]


def enough_primes(N):
    """
    Return a list of primes 2, 3, 5, ... long enough that we never run out
    of primes while enumerating candidates for a search bound N.

    The key fact used here is that, among all non-increasing exponent
    sequences with product <= N, the sequence that uses the *most*
    distinct primes is the all-ones sequence 2*3*5*7*...*p_k (the
    primorial). Any sequence with a larger exponent on an earlier prime
    only uses up the budget faster and therefore needs fewer, not more,
    distinct primes to stay under N. So it suffices to keep adding primes
    until their running product (the primorial) first exceeds N.
    """
    limit = 100
    while True:
        primes = sieve_primes(limit)
        product = 1
        for count, p in enumerate(primes, start=1):
            product *= p
            if product > N:
                return primes[:count]
        # Our sieve limit did not contain enough primes yet; grow it and retry.
        limit *= 2


def floor_log(base, value):
    """
    Return floor(log_base(value)) computed exactly, for integers base > 1
    and value >= 1.

    We start from a fast floating-point estimate and then correct it with
    exact integer power comparisons, so the result is always exactly
    right (never off by one due to floating-point rounding) while still
    being fast even when 'value' is a very large integer.
    """
    est = int(math.log(value) / math.log(base))
    while base ** est > value:
        est -= 1
    while base ** (est + 1) <= value:
        est += 1
    return est


def search_superabundant(N):
    """
    Return the sorted list of (n, sigma(n)) pairs for every superabundant
    number n <= N.

    We first enumerate every candidate n <= N of the form
    2^a1 * 3^a2 * 5^a3 * ... with a1 >= a2 >= a3 >= ... >= 1 (Theorem 1
    guarantees every superabundant number is among these candidates).
    While doing so, we drop candidates whose most recently placed
    (i.e. largest) prime exponent exceeds 1, unless the candidate is 4 or
    36 (Theorem 3 guarantees no other such candidate can be superabundant,
    so this is a safe reduction, not just a heuristic).

    Finally we sort the surviving candidates by n and scan them in order,
    keeping only the ones whose ratio sigma(n)/n strictly improves on the
    best ratio seen so far. The ratio comparison is done by cross
    multiplication of exact integers, never floating point, so there is
    no risk of rounding errors even for very large n.
    """
    primes = enough_primes(N)
    candidates = [(1, 1)]  # n = 1 is trivially superabundant

    def recurse(idx, max_exp, n, sigma):
        if idx >= len(primes):
            return
        p = primes[idx]
        power = 1  # will hold p**exp, built up incrementally
        term = 1   # will hold 1 + p + p**2 + ... + p**exp
        for exp in range(1, max_exp + 1):
            power *= p
            new_n = n * power
            if new_n > N:
                break  # increasing exp only makes new_n larger; stop here
            term += power
            new_sigma = sigma * term  # sigma is multiplicative, p is a new prime factor

            # Theorem 3 filter: this candidate's largest prime factor is p,
            # with exponent 'exp'. Keep it only if exp == 1, or if it is
            # one of the two known exceptions (4 and 36).
            if exp == 1 or new_n in (4, 36):
                candidates.append((new_n, new_sigma))

            # Continue the search regardless of the filter above: a later,
            # smaller exponent on a further prime may still yield a
            # genuinely new superabundant number.
            recurse(idx + 1, exp, new_n, new_sigma)

    initial_max_exp = N.bit_length()  # a safe upper bound: 2**k <= N needs k <= log2(N)
    recurse(0, initial_max_exp, 1, 1)

    candidates.sort()

    result = []
    best_n, best_sigma = None, None
    for n, sigma in candidates:
        if best_n is None or sigma * best_n > best_sigma * n:
            result.append((n, sigma))
            best_n, best_sigma = n, sigma
    return result


def factorize_with_known_primes(n, primes):
    """
    Factor n using trial division against a known, sufficient list of
    primes. This is only ever called on the small final list of
    superabundant numbers, so a simple trial division is fast enough
    even though n itself may be very large.
    """
    factors = []
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            exp = 0
            while n % p == 0:
                n //= p
                exp += 1
            factors.append((p, exp))
    if n > 1:
        factors.append((n, 1))
    return factors


def format_factorization(n, primes):
    """Format n as '2^a * 3^b * ...', matching the style used in the search log."""
    if n == 1:
        return "1"
    factors = factorize_with_known_primes(n, primes)
    parts = [f"{p}^{e}" if e > 1 else f"{p}" for p, e in factors]
    return " * ".join(parts)


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    if p < 1:
        raise ValueError("p must be a positive integer.")

    N = 10 ** p
    print(f"Searching for superabundant numbers up to 10^{p} = {N}")
    print()

    results = search_superabundant(N)
    primes = enough_primes(N)  # reused here only for printing factorizations

    for n, sigma in results:
        ratio = sigma / n
        factorization = format_factorization(n, primes)
        print(f"{n} = {factorization}   (sigma/n = {sigma}/{n} = {ratio})")

    print()
    print(f"Found {len(results)} superabundant numbers up to 10^{p}.")


if __name__ == "__main__":
    main()
