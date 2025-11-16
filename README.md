# p-adic-experiments
Some python scripts for doing fun things with p-adic numbers.

## What is a p-adic number?

The $p$-adic number field is a field extension of the rational numbers in the same way that the real numbers is a field extension of the rational numbers. The $p$-adic number field has its own system for writing numbers that is similar but different from that of the real numbers. The real numbers use the well-known base-10 expansion, or more generally the base-$b$ expansion, whereas $p$-adics use the $p$-adic expansion. Unlike how $b$ may be composite in the base-$b$ expansion for real numbers, it is necessary for $p$ to be a prime in the $p$-adic number field, and every $p$ creates a different number system (although the subsets of rationals isomorphic to each other).

Each digit of an integer expressed in base-$b$ gives you a little bit of information about the integer. How much does changing a single digit affect the value of the number? In base-$b$, if you change a digit that is further to the left, then the number changes "more" than if you change a digit that is further to the right. So in a sense, the digits to the left have a bigger effect on the numbers, and are more important for you to get a feel of the number. Note how this implicitly assumes the Euclidean metric $|a-b|$ to measure distances between numbers which is how we recognize a "significant" difference.

Each digit of the $p$-adic expansion of a number gives you *modular* information about the number, rather than *magnitude* information about the number. For example, the first digit tells you the equivalence class of the number modulo $p$. This is like how in base-2, the ones place digit tells you whether the number is even or odd. Together, the first two significant digits tell you the equivalence class of the number modulo $p^2$. This is analogous to how the last two digits of a non-negative integer in base-10 tell you the residue of the integer modulo 100.

From the $p$-adic perspective, we say two integers are more similar if they are equivalent modulo larger powers of $p$. Thus, the digits to the right have more significance, because they determine the congruence class of the number for more powers of $p$. For example, if you change the one's place of a 2-adic number, then the residue of that number modulo every non-negative power of two will change, but if you change a digit further to the left, then numbers residue modulo the powers of 2 below where the change happens will be unaffected.

For *non-negative integers*, the $p$-adic expansion of a number is exactly the same as the base-$p$ expansion of the number. For example, the 2-adic expansion and base-2 expansion of "67" are both `1000011.` (with infinite zeroes to the left and right).

For *negative integers*, it gets a little stranger. The $p$-adic expansion for -1 is an infinite sequence of $p-1$. For example, the 2-adic expansion of -1 is ...111111 (growing to the left) because $-1 \equiv 1 \pmod{2}$, $-1 \equiv 3 = 11_2 \pmod{4}$, etc. If you've worked with two's complement, you might recognize this. Indeed, two's complement representation is just a truncation of the 2-adic expansion: a 32-bit signed integer simply stores the residue of the integer modulo $2^{32}$, which is why algebraic operations like addition and multiplication just work (these operations respect modular congruences, since taking the modular congruence class is a [ring homomorphism](https://en.wikipedia.org/wiki/Ring_homomorphism)).

Earlier I claimed that the p-adics are a field extension of the rational numbers. But what is the equivalence class of a fraction modula a prime power? How does it make sense to ask about the residue of $2 / 3 \pmod {5}$ That unfortunately involves a little bit more complex algebra, but may make intuitive sense if you know about [Finite Fields](https://en.wikipedia.org/wiki/Finite_field).

Just like how the rational numbers make up a very small subset of all real numbers, they also make up a very small subset of all $p$-adic numbers. The $p$-adic numbers have the same cardinality as the real numbers. So when looking at the whole set, there are uncomputable $p$-adics just like how there are uncomputable real numbers.

## padic\_digitwise\_op.py

This script computes $p$-adic digitwise sum (although the code supports any digitwise operation $\textrm{op}:(\mathbb{Z}/p\mathbb{Z})\times(\mathbb{Z}/p\mathbb{Z})\longrightarrow(\mathbb{Z}/p\mathbb{Z})$ of residues modulo $p$). When $p$ is 2, the digitwise sum corresponds exactly to the bitwise XOR operation. A programmer may be familiar with the fact that XOR makes sense on non-negative integers (and some might know that it's well defined on all integers) but probably few know that *it's actually possible to compute the XOR on rational numbers too!*

Over the integers, XOR is the *unique* operation that satisfies the following properties:

 * It's an [Abelian group](https://en.wikipedia.org/wiki/Abelian_group) in which every element is an [involution](https://en.wikipedia.org/wiki/Involution_(mathematics)).
   * This means the XOR operation is associative, has an identity, is commutative, and every object is it's own inverse.
 * It respects congruence mod $2^k$ for all $k \in \mathbb{N}$.
   * This means that if $a \equiv a^\prime \pmod{2^k}$ and $b \equiv b^\prime \pmod{2^k}$, then $a \oplus b \equiv a^\prime \oplus b^\prime \pmod{2^k}$. Or in other words, it is sufficient to know $a \pmod{2^k}$ and $b \pmod{2^k}$ in order to determine $a \oplus b \pmod{2^k}$.

So in a sense, the above properties fully characterize XOR over the integers. And indeed, these properties trivially continue to hold when extending XOR to the 2-adics, which lends credibility to this generalization.

We use $p$-adic numbers rather than base-$p$ expansion ([as is done here](https://codegolf.stackexchange.com/q/147993)) to compute digitwise operations because some numbers may have multiple representations in base-$p$, and there is no good way to choose a base-$p$ representation that preserves the aforementioned XOR properties. For example, you may be familiar with how 0.9999... = 1.0000... = 1 in base-10. While, for example, you may choose to exclude base-2 expansions with finitely many one bits, and this does produce a unique representation for all rational numbers, it also breaks closure under XOR, because 0.1111... XOR 0.1111... = 0.0000... which falls into the excluded set. There isn't really a clean way to resolve this issue because solving these kind of issues is not in the domain of the base-$p$ expansion, which is more concerned with magnitudes.

### Examples

```
$ python3 padic_digitwise_op.py 1 4 2 3 2
A is 1/4; as a 2-adic, that's (0).01
B is 2/3; as a 2-adic, that's (01)10
The digit-wise sum A ⨁ B is (01)10.01
As a ratio, A ⨁ B is 11/12
```

The parentheses show the periodic repeating part that goes to the left. So above, $\frac{1}{4}$ gets represented as `(0).01` because it is equal to $2^{-2}$.

#### Associativity

Example: $\frac{2}{3} \oplus (\frac{4}{5} \oplus \frac{6}{7}) = (\frac{2}{3} \oplus \frac{4}{5}) \oplus \frac{6}{7}$
```
$ python3 padic_digitwise_op.py 4 5 6 7 2
A is 4/5; as a 2-adic, that's (0110)100
B is 6/7; as a 2-adic, that's (010)10
The digit-wise sum A ⨁ B is (010000101111)110
As a ratio, A ⨁ B is 254/65

$ python3 padic_digitwise_op.py 2 3 254 65 2
A is 2/3; as a 2-adic, that's (01)10
B is 254/65; as a 2-adic, that's (010000101111)110
The digit-wise sum A ⨁ B is (111010000101)000
As a ratio, A ⨁ B is -472/65

$ python3 padic_digitwise_op.py 2 3 4 5 2
A is 2/3; as a 2-adic, that's (01)10
B is 4/5; as a 2-adic, that's (0110)100
The digit-wise sum A ⨁ B is (1100)010
As a ratio, A ⨁ B is -22/5

$ python3 padic_digitwise_op.py -22 5 6 7 2
A is -22/5; as a 2-adic, that's (1100)010
B is 6/7; as a 2-adic, that's (010)10
The digit-wise sum A ⨁ B is (111010000101)000
As a ratio, A ⨁ B is -472/65
```

#### Identity

Example: $\frac{6}{7} \oplus 0 = 0 \oplus \frac{6}{7} = \frac{6}{7}$

```
$ python3 padic_digitwise_op.py 6 7 0 1 2
A is 6/7; as a 2-adic, that's (010)10
B is 0; as a 2-adic, that's (0)
The digit-wise sum A ⨁ B is (010)10
As a ratio, A ⨁ B is 6/7

$ python3 padic_digitwise_op.py 0 1 6 7 2
A is 0; as a 2-adic, that's (0)
B is 6/7; as a 2-adic, that's (010)10
The digit-wise sum A ⨁ B is (010)10
As a ratio, A ⨁ B is 6/7
```

#### Self-inverse

Example: $\frac{6}{7} \oplus \frac{6}{7} = 0$

```
$ python3 padic_digitwise_op.py 6 7 6 7 2
A is 6/7; as a 2-adic, that's (010)10
B is 6/7; as a 2-adic, that's (010)10
The digit-wise sum A ⨁ B is (0)
As a ratio, A ⨁ B is 0
```

#### Commutativity

Example: $\frac{6}{7} \oplus (-\frac{4}{3}) = (-\frac{4}{3}) \oplus \frac{6}{7}$

```
$ python3 padic_digitwise_op.py 6 7 -4 3 2
A is 6/7; as a 2-adic, that's (010)10
B is -4/3; as a 2-adic, that's (10)0
The digit-wise sum A ⨁ B is (000111)10
As a ratio, A ⨁ B is 14/9

$ python3 padic_digitwise_op.py -4 3 6 7 2
A is -4/3; as a 2-adic, that's (10)0
B is 6/7; as a 2-adic, that's (010)10
The digit-wise sum A ⨁ B is (000111)10
As a ratio, A ⨁ B is 14/9
```

