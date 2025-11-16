#!/bin/python3
#
# A rational number is represented as the tuple (p, q) with
# - n being the integer numerator and
# - d being the positive integer denominator.
#
# A p-adic number is represented as the tuple (seq, rpt, shr) with
# - seq being a sequence of digits in [0, p),
# - rpt being the "repeat index" (which index to cycle back to after hitting
#   the end of the sequence), and
# - shr being non-negative right-shift of seq[0] (essentially how many digits
#   passed the dot symbol are stored, which always be zero in simplified form
#   if seq[0] is zero).
# - Note that the period of repetition is len(seq) - rpt, this is very useful.

import string
import random
import math

# symbols to use to represent digits
syms = string.digits + string.ascii_lowercase

# Miller-Rabin primality test
def is_probably_prime(n, k=1000):
    if n <= 3:
        return n > 1
    elif n % 2 == 0:
        return False
    else:
        d, s = n - 1, 0
        # compute s and odd d such that 2^s * d = n - 1
        while d % 2 == 0:
            d //= 2
            s += 1
        for _ in range(k):
            a = random.randint(2, n - 2)
            x = pow(a, d, n)
            # if a^d is 1 or -1 (mod n), it's useless to us
            if x == 1 or x == n - 1:
                continue
            for _ in range(0, s - 1):
                x = pow(x, 2, n)
                # x^2 = 1 (mod n) implies x = 1 or x = -1 (mod n) for all x
                # when n is prime, but we know that x was not 1 or -1 (mod n)
                # before squaring it.
                if x == 1:
                    return False
                # if x^2 = -1, we can't learn much because after this we just
                # get (-1)^2 = 1 which holds regardless of primality.
                if x == n - 1:
                    break
            # at this point, we took our original a^d and squared it s-1 times
            # to get a^(d * 2^(s-1)). squaring it once more gets a^(d * 2^s)
            # which is a^(p-1).
            # Fermat's little theorem: x^(p-1) = 1 (mod n) for all x when n is
            # prime.
            if pow(x, 2, n) != 1:
                return False
    return True

def padic_from_rational(rat, p=2):
    n, d = rat
    # we want n / d = padic, i.e. n - padic * d = 0. this is acheived by
    # solving the equation for the padic mod p^k, one power of p at a time.
    seq = []
    seen = {}
    shr = 0
    # first let's make gcd(d, p) = 1
    while d % p == 0:
        d //= p
        if n % p == 0:
            # this simplification step helps us avoid unnecessary trailing
            # zeroes
            n //= p
        else:
            # we're removing a p-factor in d but not in n, so we need to shift
            # the resulting padic to have the same effect of dividing it by p.
            shr += 1
    di = pow(d, -1, p)
    # induction on k: want to preserve n - curr_padic * d = 0 (mod p^k)
    # - we'll just store (n - curr_padic * d) / (p^k) as our state in "n"
    # - we want to add a digit to curr_padic to solve that initial congruence
    #   mod p^(k+1).
    #   - we can do this by adding a digit to seq in the p^k place, namely our
    #     current "n" state times the modular multiplicative inverse of "d".
    # - if we come across an "n" state we have seen before, we know the future
    #   outcomes will be the same, so we found our repeat index.
    while n not in seen:
        # store the repeat index len(seq) for the state value n for future use
        seen[n] = len(seq)
        # i is the digit we need to add to solve the congruence
        i = (n * di) % p
        seq.append(i)
        # update the state n with the new digit in the padic, making it now
        # congruent to zero mod p.
        n -= d * i
        # divide the state by p to setup for next iteration
        n //= p
    return seq, seen[n], shr

def rational_from_padic(padic, p=2):
    seq, rpt, shr = padic
    # "k" keeps track of the power of p for the current digit
    a, c, k = 0, 0, 1
    # accumulate the non-repeating digits into an integer "a"
    for i in range(rpt):
        a += seq[i] * k
        k *= p
    # for the repeating digits, we can use the geometric series formula:
    #   c + c * r + c * r^2 ... = c / (1 - r).
    # if the repeating part repeats with period t, and has a value of c in the
    # first iteration, then it's equal to
    #   c + c * p^t + c * p^2t ...
    # so we apply the geometric series formula with c = c and r = p^t to get
    #   c / (1 - p^t).
    for i in range(rpt, len(seq)):
        c += seq[i] * k
        k *= p
    # we can apply shr by just multiplying by p^-shr. finally, the number is
    # (a + c / (1 - p^t)) * p^-shr which expands into a fraction as
    # (a * (p^t - 1) - c) / ((p^t - 1) * p^shr)
    d = pow(p, len(seq) - rpt) - 1 # d = p^t - 1
    return a * d - c, d * pow(p, shr)

def simplify_rational(rat):
    n, d = rat
    # make denominator positive
    if d < 0:
        n, d = -n, -d
    # make numerator and denominator coprime
    g = math.gcd(n, d)
    return n // g, d // g

def simplify_padic(padic):
    seq, rpt, shr = padic
    # t is the number of trailing zeroes past the dot symbol
    t = next((i for (i, d) in enumerate(seq[:shr]) if d != 0), min(len(seq), shr))
    if t == len(seq): # all zeroes
        return [0], 0, 0
    seq, rpt, shr = seq[t:], rpt - t, shr - t
    # shift the repeated part as far to the right as possible by merging in the
    # non-repeated digits if they match up.
    while rpt > 0 and seq[rpt - 1] == seq[-1]:
        rpt -= 1
        seq.pop()
    # check if the repeated part actually has a smaller period
    for i in range(1, (len(seq) - rpt) // 2 + 1):
        # if i divides the period...
        if (len(seq) - rpt) % i == 0:
            # is the repeated part periodic with period i?
            if all(seq[j] == seq[j-i] for j in range(rpt + i, len(seq))):
                # if so, then truncate it!
                seq = seq[:rpt+i]
                # since we iterate the factors of the period in ascending
                # order, all periods less than i have already been checked, so
                # we can eagerly break the loop.
                break
    return seq, rpt, shr

def str_from_rational(rat):
    n, d = rat
    if d == 1:
        return f"{n}"
    else:
        return f"{n}/{d}"

def str_from_padic(padic):
    seq, rpt, shr = padic
    if max(seq) >= len(syms):
        return "{not enough symbols to represent this number}"
    # in the repeating part, we want to print exactly len(seq) - rpt digits.
    # we want the last digit to be max(shr, rpt) because we don't want to print
    # the non-repeating digits or the digits right of the dot symbol here.
    s = "("
    for i in range((max(shr - rpt, 0) + len(seq)) - 1, max(shr, rpt) - 1, -1):
        s += syms[seq[rpt + ((i - rpt) % (len(seq) - rpt))]]
    s += ")"
    # we can print the non-repeating digits left of the dot symbol, if any
    for i in range(max(shr, rpt) - 1, shr - 1, -1):
        s += syms[seq[i]]
    if shr:
        # we can print the symbols to the right of the dot symbol.
        # these might not all be present in seq if seq is simplified.
        # for example: 1/16 with p=2 is stored as seq=[1, 0], rpt=1, shr=4.
        s += "."
        for i in range(shr - 1, min(shr, rpt) - 1, -1):
            s += syms[seq[rpt + ((i - rpt) % (len(seq) - rpt))]]
        for i in range(min(shr, rpt) - 1, -1, -1):
            s += syms[seq[i]]
    return s

def padic_digitwise_op(op, padic_a, padic_b):
    seen = {}
    seq = []
    seq_a, rpt_a, shr_a = padic_a
    seq_b, rpt_b, shr_b = padic_b
    # idx_a and idx_b are indices in seq_a and seq_b. we need to start with the
    # logical digit position corresponding to the larger shr value. for example
    # if shr_a = 5 and shr_b = 3, we start with idx_a = 0 and idx_b = -2.
    idx_a, idx_b = min(0, shr_a - shr_b), min(0, shr_b - shr_a)
    # the repetition state that we need to track is the pair (idx_a, idx_b). if
    # these indices repeat, then the future outcomes will be the same.
    while (idx_a, idx_b) not in seen:
        seen[idx_a, idx_b] = len(seq)
        # get the current logical digit from both sequences. a negative index
        # means that other sequence had a larger shr value, we can just use
        # zero.
        dig_a = seq_a[idx_a] if idx_a >= 0 else 0
        dig_b = seq_b[idx_b] if idx_b >= 0 else 0
        # apply the op. the first digit in the resulting seq is in the logical
        # -max(shr_a, shr_b) position.
        seq.append(op(dig_a, dig_b))
        # advance the logical indices, repeating if necessary.
        idx_a += 1
        idx_b += 1
        if idx_a == len(seq_a):
            idx_a = rpt_a
        if idx_b == len(seq_b):
            idx_b = rpt_b
    return seq, seen[idx_a, idx_b], max(shr_a, shr_b)

sum_op = lambda p: lambda a, b: (a + b) % p

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print(f"Usage: {sys.argv[0]} <numerA> <denomA> <numerB> <denomB> [p;default=2]")
        print(f"Computes the digit-wise sum of numerA/denomA and numerB/denomB in the p-adic field")
        print(f"For p=2, this is the XOR of the two numbers")
        sys.exit(1)
    rat_a = simplify_rational((int(sys.argv[1]), int(sys.argv[2])))
    rat_b = simplify_rational((int(sys.argv[3]), int(sys.argv[4])))
    p = 2 if len(sys.argv) < 6 else int(sys.argv[5])
    assert is_probably_prime(p), "p must be prime"
    assert rat_a[1] > 0, "denomA must be positive"
    assert rat_b[1] > 0, "denomB must be positive"
    padic_a = padic_from_rational(rat_a, p=p)
    padic_b = padic_from_rational(rat_b, p=p)
    padic_dws = padic_digitwise_op(sum_op(p), padic_a, padic_b)
    padic_dws = simplify_padic(padic_dws)
    rat_dws = simplify_rational(rational_from_padic(padic_dws, p))
    print(f"A is {str_from_rational(rat_a)}; as a {p}-adic, that's {str_from_padic(padic_a)}")
    print(f"B is {str_from_rational(rat_b)}; as a {p}-adic, that's {str_from_padic(padic_b)}")
    print(f"The digit-wise sum A ⨁ B is {str_from_padic(padic_dws)}")
    print(f"As a ratio, A ⨁ B is {str_from_rational(rat_dws)}")
