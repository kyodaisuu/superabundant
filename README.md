# Superabundant Number Search

A Python program for searching for **superabundant numbers** using the structural theorems of Alaoglu and Erdős (1944).

The program exploits the special structure of superabundant numbers to search far more efficiently than testing every integer up to a given bound. It uses arbitrary-precision integers provided by Python, allowing searches beyond the range where conventional sieve-and-factor approaches become practical.

## Definition

Let σ(n) denote the sum of the positive divisors of n:

```
σ(n) = sum of d over all positive divisors d of n.
```

A positive integer n is called **superabundant** if

```
σ(m) / m < σ(n) / n
```

for every integer m satisfying

```
1 <= m < n.
```

In other words, n is superabundant if the ratio σ(n)/n reaches a new record at n.

The program searches for all superabundant numbers up to a specified bound.

## Mathematical Background

The search is based on three structural theorems proved by Leonidas Alaoglu and Paul Erdős in their 1944 paper:

> Alaoglu, L. and Erdős, P. (1944),
> "On Highly Composite and Similar Numbers",
> Transactions of the American Mathematical Society, 56, 448–469.

### Theorem 1: Non-increasing exponents

If

```
n = 2^a1 * 3^a2 * 5^a3 * ... * p_k^ak
```

is superabundant, then

```
a1 >= a2 >= a3 >= ... >= ak >= 1.
```

Furthermore, the prime factors must be consecutive primes starting from 2. Thus, a superabundant number can only have the form

```
2^a1 * 3^a2 * 5^a3 * ... * p_k^ak
```

with non-increasing exponents.

This dramatically reduces the search space compared with checking every integer.

### Theorem 2: Exponents of later primes

Let q < r be primes, and let a_q be the exponent of q. Define

```
beta = floor(a_q * log(q) / log(r)).
```

Then the exponent a_r of r in a superabundant number must satisfy

```
a_r in {beta - 1, beta, beta + 1}.
```

In this implementation, the theorem is applied with q = 2.

Once the exponent a1 of 2 is chosen, the theorem therefore restricts the exponent of every subsequent prime to at most three possible values. This reduces the branching factor of the search from potentially `a_(i-1)` possibilities to at most three.

### Theorem 3: Exponent of the largest prime

If p is the largest prime factor of a superabundant number, then the exponent of p is 1, with only two exceptions:

```
n = 4
n = 36
```

The program uses this theorem to discard many candidates during the search.

## Search Algorithm

The program proceeds in three main stages.

### 1. Enumerate structurally possible candidates

Using Theorem 1, the program recursively enumerates numbers of the form

```
2^a1 * 3^a2 * 5^a3 * ... * p_k^ak
```

where

```
a1 >= a2 >= ... >= ak >= 1.
```

The search is restricted to numbers satisfying

```
n <= N.
```

Theorem 2 is used to restrict the exponent of each newly introduced prime to at most three candidate values.

Theorem 3 is also applied during the search to eliminate candidates whose largest prime factor has an exponent greater than 1, except for 4 and 36.

These are mathematically valid reductions, not heuristic pruning: every superabundant number is guaranteed to survive these structural constraints.

### 2. Compute σ(n)

For each surviving candidate, σ(n) is computed using the multiplicativity of the divisor-sum function.

For

```
n = p1^a1 * p2^a2 * ... * pk^ak,
```

we have

```
σ(n) = (1 + p1 + ... + p1^a1)
     * (1 + p2 + ... + p2^a2)
     * ...
     * (1 + pk + ... + pk^ak).
```

The program computes these values incrementally while constructing candidates.

### 3. Identify record-holders

The surviving candidates are sorted by n.

They are then scanned in increasing order, and a candidate is reported as superabundant if

```
σ(n) / n
```

is strictly larger than the best ratio seen so far.

The comparison is performed using exact integer cross multiplication:

```
σ(n1) / n1 < σ(n2) / n2
```

is tested as

```
σ(n1) * n2 < σ(n2) * n1.
```

Thus, the identification of superabundant numbers does not depend on floating-point arithmetic or numerical rounding.

## Prime Generation

The program determines how many primes are needed for the search using the primorial

```
2 * 3 * 5 * 7 * ... * p_k.
```

Among all non-increasing exponent sequences, the sequence consisting entirely of exponents 1 uses the smallest possible exponent on each prime and therefore allows the largest possible number of distinct prime factors under a fixed bound N.

Consequently, the program generates primes until the primorial first exceeds N. This guarantees that enough primes are available to enumerate every structurally possible candidate up to N.

## Arbitrary-Precision Integers

Python integers have arbitrary precision, so the program does not have a fixed machine-integer limit on the size of N.

The command-line argument specifies an exponent p, and the search bound is

```
N = 10^p.
```

For example:

```
python superabundant.py 16
```

searches all superabundant numbers up to

```
10^16.
```

In principle, the program can accept arbitrarily large values of p, although the actual computation time and memory requirements naturally increase with the size of the search space.

## Requirements

* Python 3
* Standard library only

No external Python packages are required.

## Usage

Run the program with:

```
python superabundant.py [p]
```

where `p` is a positive integer. The program searches up to `10^p`.

For example:

```
python superabundant.py 16
```

The default value is `16`, so simply running

```
python superabundant.py
```

is equivalent to:

```
python superabundant.py 16
```

The output lists each superabundant number, its prime factorization, and its σ(n)/n ratio.

For example, the output has the following general form:

```
1 = 1   (sigma/n = 1/1 = 1.0)
2 = 2   (sigma/n = 3/2 = 1.5)
4 = 2^2   (sigma/n = 7/4 = 1.75)
...
```

At the end, the program reports the total number of superabundant numbers found.

## Example

To search up to 10^40:

```
$ python superabundant.py 40
```

The program prints:

```
Searching for superabundant numbers up to 10^40 = 10000000000000000000000000000000000000000
```

followed by the detected superabundant numbers and their factorizations.

## Why This Is Efficient

A naive approach would examine every integer

```
1, 2, 3, ..., N
```

and determine whether each one is superabundant. This quickly becomes impractical because it requires factoring a huge number of integers and computing their divisor sums.

This program instead searches only a highly restricted family of candidates:

1. **Theorem 1** restricts prime factors to consecutive primes and requires non-increasing exponents.
2. **Theorem 2** reduces the possible exponent of each new prime to at most three values.
3. **Theorem 3** eliminates candidates whose largest prime factor has an exponent greater than 1, apart from 4 and 36.
4. Exact integer arithmetic avoids expensive or unreliable high-precision floating-point comparisons.

The combination of these structural properties makes it possible to search much larger ranges than a straightforward exhaustive search based only on Theorem 1.

## Limitations

Although the algorithm is substantially more efficient than exhaustive enumeration, the search is still computationally demanding for extremely large bounds.

The number of structurally possible candidates grows as the bound increases, and the final candidate list must be sorted before the record-holders can be identified.

The implementation is therefore intended primarily as a computational exploration of superabundant numbers and the structural results of Alaoglu and Erdős, rather than as an asymptotically optimal algorithm for computing arbitrarily large superabundant numbers.

## Result

Result for calculating up to 10^100 is saved at `result-100.txt`.

## Reference

Alaoglu, L. and Erdős, P. (1944). "On Highly Composite and Similar Numbers." Transactions of the American Mathematical Society, 56, 448–469.

## License

MIT License
