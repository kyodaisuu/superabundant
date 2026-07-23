#!/usr/bin/env python3
"""
Superabundant number search using the structural theorems of
Alaoglu and Erdos (1944), "On Highly Composite and Similar Numbers",
Trans. Amer. Math. Soc. 56, 448-469.

A number n is superabundant if sigma(m)/m < sigma(n)/n for all 0 < m < n,
where sigma(n) is the sum of divisors of n.

Theorem 1 (Alaoglu-Erdos):
If

    n = 2^a1 * 3^a2 * 5^a3 * ... * p_k^ak

is superabundant, then

    a1 >= a2 >= a3 >= ... >= ak >= 1.

In particular, a superabundant number uses consecutive primes starting
from 2, with non-increasing exponents.

Theorem 2 (Alaoglu-Erdos):
Let q < r be primes, and let a_q and a_r be their exponents in a
superabundant number. Set

    beta = floor(a_q * log(q) / log(r)).

Then

    a_r in {beta - 1, beta, beta + 1}.

This implementation applies Theorem 2 for every previously selected
prime q < r, not just q = 2. The possible exponents for r are therefore
restricted to the intersection of all such three-element sets, together
with the non-increasing-exponent constraint from Theorem 1.

Theorem 3 (Alaoglu-Erdos):
If p is the largest prime factor of a superabundant number, then the
exponent of p is 1, except for n = 4 and n = 36.

Theorem 3 is used to exclude candidates from the final list during the
search. However, the recursion is continued from such candidates,
because adding a new prime factor with exponent 1 can produce a number
whose largest prime exponent satisfies Theorem 3.

The final identification of superabundant numbers is performed by
sorting all surviving candidates by n and retaining those for which
sigma(n)/n is a new record.

All comparisons of sigma(n)/n use exact integer cross multiplication.
"""

import math
import sys


def sieve_primes(limit):
    """Return the list of primes <= limit using a sieve of Eratosthenes."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime[i]:
            is_prime[
                i * i : limit + 1 : i
            ] = bytearray(len(range(i * i, limit + 1, i)))

    return [i for i, flag in enumerate(is_prime) if flag]


def enough_primes(N):
    """
    Return enough consecutive primes to enumerate every candidate <= N
    allowed by Theorem 1.

    The maximum possible number of distinct prime factors is obtained
    when every exponent is 1, giving the primorial

        2 * 3 * 5 * 7 * ...

    We therefore generate primes until the primorial first exceeds N.
    """
    limit = 100

    while True:
        primes = sieve_primes(limit)
        product = 1

        for count, p in enumerate(primes, start=1):
            product *= p

            if product > N:
                return primes[:count]

        limit *= 2


def floor_log_ratio(a, q, r):
    """
    Return floor(a * log(q) / log(r)) exactly.

    This is equivalent to finding the largest integer k such that

        r^k <= q^a.

    A floating-point estimate is used as a starting point, then corrected
    using exact integer comparisons.

    Parameters
    ----------
    a : int
        Exponent of q.
    q : int
        Smaller prime.
    r : int
        Larger prime.
    """
    if a < 0:
        raise ValueError("a must be non-negative.")

    if q <= 1 or r <= 1:
        raise ValueError("q and r must be greater than 1.")

    if q >= r:
        raise ValueError("q must be smaller than r.")

    if a == 0:
        return 0

    # Initial floating-point estimate.
    estimate = int(a * math.log(q) / math.log(r))

    # Exact correction.
    # We need the largest k satisfying r^k <= q^a.
    q_power = q ** a

    while estimate > 0 and r ** estimate > q_power:
        estimate -= 1

    while r ** (estimate + 1) <= q_power:
        estimate += 1

    return estimate


def theorem2_exponent_candidates(exponents, primes, r_index):
    """
    Return the possible exponent values for primes[r_index] according
    to Theorem 2 applied to every previously selected prime.

    If the exponent of q is a_q, Theorem 2 gives

        a_r in {
            beta - 1,
            beta,
            beta + 1
        }

    where

        beta = floor(a_q * log(q) / log(r)).

    The possible values from all q < r are intersected.

    Theorem 1 is not applied here; its constraint
    a_r <= a_(r-1) is imposed separately by the caller.
    """
    r = primes[r_index]

    possible = None

    for q_index in range(r_index):
        q = primes[q_index]
        a_q = exponents[q_index]

        beta = floor_log_ratio(a_q, q, r)

        candidates = {
            beta - 1,
            beta,
            beta + 1,
        }

        # Exponents must be positive.
        candidates = {
            exponent
            for exponent in candidates
            if exponent >= 1
        }

        if possible is None:
            possible = candidates
        else:
            possible &= candidates

        # Once the intersection is empty, no exponent is possible.
        if not possible:
            return []

    return sorted(possible)


def search_superabundant(N):
    """
    Return the sorted list of (n, sigma(n)) pairs for every superabundant
    number n <= N.

    Candidate generation uses:

    1. Theorem 1:
       Prime factors are consecutive and exponents are non-increasing.

    2. Theorem 2:
       For every new prime r, its exponent is restricted by every
       previously selected prime q < r. The resulting candidate sets
       are intersected.

    3. Theorem 3:
       A candidate whose largest prime factor has exponent greater than 1
       cannot itself be superabundant, except for 4 and 36.

    Theorem 3 does not terminate recursion from a rejected candidate,
    because a later prime factor may restore the condition that the
    largest prime has exponent 1.

    The final list is obtained by sorting candidates by n and retaining
    exactly those for which sigma(n)/n is a new record.
    """
    primes = enough_primes(N)

    # n = 1 is superabundant by definition.
    candidates = [(1, 1)]

    def recurse(idx, exponents, n, sigma):
        """
        Add primes[idx], primes[idx + 1], ... recursively.

        `exponents` contains the exponents of all previously selected
        primes.
        """
        if idx >= len(primes):
            return

        p = primes[idx]

        # Theorem 2:
        # Apply it to every previously selected prime q < p.
        candidate_exps = theorem2_exponent_candidates(
            exponents,
            primes,
            idx,
        )

        if not candidate_exps:
            return

        # Theorem 1:
        # Exponents must be non-increasing.
        max_exp = exponents[-1]

        candidate_exps = [
            exp
            for exp in candidate_exps
            if exp <= max_exp
        ]

        for exp in candidate_exps:
            power = p ** exp
            new_n = n * power

            if new_n > N:
                continue

            # sigma(p^exp) = 1 + p + ... + p^exp
            term = (power * p - 1) // (p - 1)

            # sigma is multiplicative because p is a new prime factor.
            new_sigma = sigma * term

            new_exponents = exponents + [exp]

            # Theorem 3:
            # The current largest prime p must have exponent 1,
            # except for n = 4 and n = 36.
            #
            # This determines whether the current candidate itself
            # can be superabundant. We nevertheless continue recursion
            # regardless of this test.
            if exp == 1 or new_n in (4, 36):
                candidates.append((new_n, new_sigma))

            recurse(
                idx + 1,
                new_exponents,
                new_n,
                new_sigma,
            )

    # Choose the exponent of 2 first.
    #
    # Theorem 2 requires a previously selected smaller prime q.
    # Therefore, once the exponent of 2 is chosen, every subsequent
    # prime can be constrained by Theorem 2.
    initial_max_exp = N.bit_length()

    for first_exp in range(1, initial_max_exp + 1):
        n = 2 ** first_exp

        if n > N:
            break

        sigma = 2 ** (first_exp + 1) - 1

        # For n = 2^a, Theorem 3 says that only a = 1 is superabundant,
        # except for n = 4.
        if first_exp == 1 or n in (4, 36):
            candidates.append((n, sigma))

        recurse(
            idx=1,
            exponents=[first_exp],
            n=n,
            sigma=sigma,
        )

    # Multiple recursion paths should not normally produce the same n,
    # but using a set here makes the final stage robust.
    candidates = sorted(set(candidates))

    result = []

    best_n = None
    best_sigma = None

    for n, sigma in candidates:
        if (
            best_n is None
            or sigma * best_n > best_sigma * n
        ):
            result.append((n, sigma))
            best_n = n
            best_sigma = sigma

    return result


def factorize_with_known_primes(n, primes):
    """
    Factor n using trial division against a known sufficient list of primes.

    This is called only for the final list of superabundant numbers.
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
    """Format n as '2^a * 3^b * ...'."""
    if n == 1:
        return "1"

    factors = factorize_with_known_primes(n, primes)

    parts = [
        f"{p}^{e}" if e > 1 else f"{p}"
        for p, e in factors
    ]

    return " * ".join(parts)


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 16

    if p < 1:
        raise ValueError("p must be a positive integer.")

    N = 10 ** p

    print(
        f"Searching for superabundant numbers up to "
        f"10^{p} = {N}"
    )
    print()

    results = search_superabundant(N)
    primes = enough_primes(N)

    for n, sigma in results:
        ratio = sigma / n
        factorization = format_factorization(n, primes)

        print(
            f"{n} = {factorization}   "
            f"(sigma/n = {sigma}/{n} = {ratio})"
        )

    print()
    print(
        f"Found {len(results)} superabundant numbers "
        f"up to 10^{p}."
    )


if __name__ == "__main__":
    main()
