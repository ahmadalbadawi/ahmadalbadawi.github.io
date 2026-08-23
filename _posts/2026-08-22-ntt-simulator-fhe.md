---
title: 'A Step-by-Step Simulator for the NTT in FHE'
date: 2026-08-22
permalink: /posts/2026/08/ntt-simulator-fhe/
mathjax: true
description: 'An interactive simulator that replays the real FHE NTT one butterfly at a time: recorded traces and the six places where the production code differs from the textbook.'
tags:
  - FHE
  - NTT
  - Number Theoretic Transform
  - OpenFHE
  - Polynomial Multiplication
  - Negacyclic Convolution
  - Homomorphic Encryption
  - Interactive Simulator
---

**[Ahmad Al Badawi](https://ahmadalbadawi.com/)**

***

In most FHE schemes, the dominant mathematical operation is polynomial multiplication in a ring. To do these multiplications efficiently, we use a mathematical trick known as the Number Theoretic Transform (NTT). If you know DFT and FFT, NTT should come naturally to you once you understand the field it works in. Most FHE libraries include a variant of NTT implementation. Taking [OpenFHE](https://github.com/openfheorg/openfhe-development/tree/v1.5.1) as an example, the code that performs NTT is two functions in one header file: [`ForwardTransformToBitReverseInPlace`](https://github.com/openfheorg/openfhe-development/blob/v1.5.1/src/core/include/math/hal/intnat/transformnat-impl.h#L302-L374) for the forward transform and [`InverseTransformFromBitReverseInPlace`](https://github.com/openfheorg/openfhe-development/blob/v1.5.1/src/core/include/math/hal/intnat/transformnat-impl.h#L511-L625) for the inverse. Neither function looks like the algorithm in a textbook.

I built a simulator that runs that code one butterfly at a time and shows what each line does: the **[NTT / INTT Explorer](https://visuallearn.github.io/fhe-ntt-sim/)**. This article introduces the simulator. It names the places where the production code and the textbook disagree.

***Disclaimer**: This article is written for readers who know what a polynomial is and what modular arithmetic is. You do not need prior experience with the NTT or with the internals of OpenFHE.*

***

**Contents**
* ToC
{:toc}


***

## What the Simulator Is

The simulator is a static web page. The browser contains no transform code at all. Every value on the screen comes from an instrumented build of OpenFHE 1.5.1, and the page replays those recordings.

It covers ring dimensions 4, 8, 16, and 32, with moduli of 5 to 10 bits. It uses one modulus, not a residue number system chain. The numbers are small enough to check with a pen.

Each step shows five things at once:

- a dataflow diagram of the whole transform, with the current butterfly marked
- the array in memory, and the two slots that the butterfly rewrites
- the twiddle factor table, with the entry in use highlighted
- the arithmetic of this one step, with the actual numbers in it
- the OpenFHE source code, with the cursor on the line that runs

There is a twelve-step guided tour for a first visit. A second view shows where the evaluation points come from. A third view compares polynomial multiplication with and without the transform.

***

## Why the Source Code Is Hard to Read

A textbook forward NTT is about fifteen lines. The one in OpenFHE is 73 lines. The extra length is not accidental. Six engineering decisions explain almost all of it, and each one stays invisible until somebody points it out.

The examples below all use $N = 8$ and modulus $q = 17$, where $\psi = 3$. The simulator uses the same configuration by default. Here, $N$ is the ring dimension, $q$ is the modulus, and $\psi$ is a primitive $2N$-th root of unity modulo $q$.

### 1. The Twiddle Table Is Not in Order

The transform needs the powers of $\psi$ many times. OpenFHE calculates them once and keeps them in a table. The table is not in natural order, but rather in a bit-reversed order:

| slot | binary | reversed | power | value |
|---|---|---|---|---|
| 0 | 000 | 000 | $\psi^0$ | 1 |
| 1 | 001 | 100 | $\psi^4$ | 13 |
| 2 | 010 | 010 | $\psi^2$ | 9 |
| 3 | 011 | 110 | $\psi^6$ | 15 |
| 4 | 100 | 001 | $\psi^1$ | 3 |
| 5 | 101 | 101 | $\psi^5$ | 5 |
| 6 | 110 | 011 | $\psi^3$ | 10 |
| 7 | 111 | 111 | $\psi^7$ | 11 |

Slot $j$ holds $\psi$ raised to the bit-reversal of $j$. For example, slot 1 is `001` in binary. The reversal is `100`, which is 4. So slot 1 holds $\psi^4 = 13$.

This order looks strange until you read the loops. The butterflies read the table in sequence: slot 1, then slot 2, then slot 3. The bit-reversed order makes that simple walk give the correct factor every time. The inner loop needs no index arithmetic at all.

### 2. The Negacyclic Twist Is Inside the Same Table

The polynomials live in a ring where $X^N = -1$. A product that goes beyond degree $N$ wraps around to a lower degree and changes sign.

A textbook handles that with a separate pass. It multiplies coefficient $i$ by $\psi^i$ before the transform, and by $\psi^{-i}$ after the inverse. OpenFHE does not have that pass. Its table holds powers of a $2N$-th root of unity, not an $N$-th root. The same butterflies then produce the twist with no extra work.

The result is a clean statement of what the forward transform calculates:

$$ \text{output}[p] = a\left(\psi^{\,2\,\text{brev}(p)+1}\right) \bmod q $$

Here, $\text{brev}(p)$ is the bit-reversal of $p$. The transform evaluates the polynomial at the **odd** powers of $\psi$. One table does two jobs.

### 3. The Output Stays Out of Order

When the last stage finishes, the array holds the $N$ values. The order is not the one you expect. Slot 0 holds the value at $\psi^1$. Slot 1 holds the value at $\psi^9$. Slot 2 holds the value at $\psi^5$.

OpenFHE does not sort this output, because only two operations read it. The first is the pointwise multiplication, which uses one product for each slot, so the order does not matter there. The second is the inverse transform, and it is written to expect exactly this order. A sort costs real time and gives nothing.

This is a common source of confusion when you read FHE code. The value form of a polynomial has a permuted order, and every later operation uses the same permutation. If you read my article on [SIMD packing in BGV/BFV](https://ahmadalbadawi.com/posts/2025/12/simd-packing-bgv-bfv-fhe/), you saw the same idea: a polynomial has a coefficient form and a value form, and each form is good at different operations.

### 4. One Stage Is a Separate Loop

The forward transform writes its last stage as its own loop, outside the main loop. The inverse transform does the same with its first stage. Engineers call this peeling. The general case then has no stride of one, and the compiler optimizes it better.

The butterflies are identical. Only the code is separate. When you reach those steps, the simulator moves the source cursor to that other loop. It also shades the affected column of the diagram and adds a note that says why.

### 5. The Division by $N$ Hides Inside a Twiddle

An inverse transform ends with a division by $N$. OpenFHE does most of that division inside the last twiddle:

```
omega1Inv = psi^-(N/2) * N^-1
          = 4 * 15  mod 17
          = 9
```

That twiddle divides the upper half of the array by $N$ as a side effect of its normal butterfly. A separate pass then multiplies only the lower half by $N^{-1}$. The saving is $N/2$ multiplications. [OpenFHE issue 872](https://github.com/openfheorg/openfhe-development/issues/872) records the reason.

The cost is readability. The last twiddle is not a power of $\psi^{-1}$, and that looks like a mistake to a reader who does not know the trick. The simulator flags that step and shows where the 9 comes from.

### 6. The Modular Multiplication Is Not a Modulo Operation

The inner loop reads:

```cpp
omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
```

Many butterflies use the same twiddle, so OpenFHE calculates `floor(w * 2^64 / q)` once. Each multiplication then needs one 64-bit multiply-high, two low multiplies, one subtraction, and one conditional addition. It needs no division, and that is the part that matters. This is the method of Shoup from NTL, analyzed in Harvey's paper on faster NTT arithmetic ([Harvey, 2014](#ref-harvey)). In other places, OpenFHE uses Barrett reduction.

These optimization tricks change the cost, but the end result stays the same.

***

## One Butterfly, With Real Numbers

Everything above is structure. Here is a single step, taken from the recorded trace of the ramp input at $N = 8$ and $q = 17$.

The array holds `[1, 15, 12, 9, 16, 4, 9, 14]` as stage 2 begins. Stage 2 uses a stride of 2 and reads table slot 2, so the twiddle is $\psi^2 = 9$. Slot 2 is its own bit-reversal, and that is why the slot number and the exponent agree here. They do not agree at slot 1, as the table above shows. The first butterfly of the stage works on slots 0 and 2:

```
u = X[0] = 1        v = X[2] = 12        omega = 9

v * omega  = 12 * 9  = 108  = 6   mod 17
X[0] <- u + v*omega  = 1 + 6  = 7   mod 17
X[2] <- u - v*omega  = 1 - 6  = -5  = 12  mod 17
```

One multiplication, one addition, one subtraction, and two slots rewritten in place. That is the whole operation. A transform at $N = 8$ performs 12 of them, in three stages of four.

The simulator shows this step in all five panels at once: the dataflow diagram, the array, the table row in use, the arithmetic, and the source line. The arithmetic above is what the working panel prints.

***

## The Parameters Are Not Free Either

A $2N$-th root of unity exists only when $q \equiv 1 \pmod{2N}$. OpenFHE selects the modulus with `LastPrime<NativeInteger>(bits, 2N)`: the largest prime that has exactly that many bits and satisfies the congruence.

The constraint removes more combinations than you expect. The simulator covers 28 pairs of ring dimension and modulus size. Only 18 of them are possible:

| bits | $N=4$ | $N=8$ | $N=16$ | $N=32$ |
|---|---|---|---|---|
| 4 | none | none | none | none |
| 5 | 17 | 17 | none | none |
| 6 | 41 | none | none | none |
| 7 | 113 | 113 | 97 | none |
| 8 | 241 | 241 | 193 | 193 |
| 9 | 457 | 449 | 449 | 449 |
| 10 | 1009 | 1009 | 929 | 769 |

A 4-bit modulus is impossible for every ring dimension here. $N = 32$ needs at least 8 bits. The simulator shows this grid and disables the cells that cannot exist, with OpenFHE's own error message on each one.

***

## What I Simplified, and Why

The working panel shows plain `(a * omega) mod q`. It does not animate the Shoup multiplication from decision 6 above.

That was a deliberate choice, and it is the one I am least sure about. At $q = 17$, the precomputed constants are nineteen or twenty digits long. The values they multiply have two digits. Those constants hide the transform instead of explaining it. The real constants are in the data, and a note on every multiplication states what OpenFHE actually does. A reader who wants the reduction technique gets the constants and the citation. A reader who wants the transform does not have to read twenty-digit numbers first.

The second simplification is one modulus, not a chain. Real OpenFHE ciphertexts use a residue number system and run this same transform once for each modulus in the chain. Here there is exactly one modulus, so no Chinese Remainder Theorem arithmetic hides the butterflies.

***

## One Number That Is Easy to Misread

The guided tour counts the cost of one polynomial multiplication through the transform. The count is three transforms, plus $N$ pointwise products, plus $N/2$ multiplications by $N^{-1}$ at the end of the inverse. At $N = 8$, that is $36 + 8 + 4 = 48$ multiplications, against 64 for the direct method. That is not impressive.

The textbook figure for this comparison is 44, because the convention counts the three transforms and the pointwise products and stops there. The simulator animates the $N^{-1}$ steps, so its own total counts them.

It is also not the real comparison. A ciphertext stays in the value form. Additions and multiplications there need one operation for each point. A program converts once, does many operations, and converts back only when it needs the coefficients. OpenFHE keeps a ciphertext in the value form for exactly this reason. The program pays the conversion cost once for a large part of a calculation, not once per multiplication.

At the ring dimensions real schemes use, the direct method is not an option at all.

***

## How to Read It

Open the **[guided tour](https://visuallearn.github.io/fhe-ntt-sim/#/tour?t=0)** first. It works through one configuration, $N = 8$ with $q = 17$, in twelve steps.

Then open the Transform view. Select a ring dimension and a modulus size. Select an input polynomial. Move with the left and right arrow keys. The up and down arrows move one whole stage. The space bar plays.

Every view has a deep link that carries its state. A specific butterfly at a specific configuration has an address that you can send to somebody.

***

## Reproducing the Data

Everything rebuilds from a pinned OpenFHE release with one script. The build needs g++, CMake, git, and python3. It needs no Node and no npm. The page is plain ES modules, and a server sends the files without changes.

```sh
git clone https://github.com/visuallearn/fhe-ntt-sim && cd fhe-ntt-sim
./tools/06_all.sh     # fetch, patch, build, generate, and run every test
```

The reference experiments build separately, against an installed OpenFHE, in the same way as any external project:

```sh
./tools/12_install_openfhe.sh
./tools/14_build_gt_exp.sh
./openfhe-gt-exp/build/gt_experiment 8 5
```

The last command prints one experiment:

- the parameters, with the OpenFHE call that produced each one
- both twiddle tables, with their Shoup constants
- the twiddles that each stage reads
- the evaluation points
- every input and output vector, compared against an independent oracle

If you find a value that looks wrong, the data is checkable and the pipeline is reproducible. Please tell me on my [contact page](/contact/), or open an issue on the [repository](https://github.com/visuallearn/fhe-ntt-sim).

***

## Key Takeaways

*   **The production NTT differs from the textbook in six specific decisions:** a bit-reversed twiddle table, the negacyclic twist folded into that table, output that stays permuted, one peeled stage, a division by $N$ hidden inside a twiddle, and Shoup modular multiplication. Each one is a small optimization. Together, they turn fifteen textbook lines into 73.
*   **The simulator replays real recordings.** The browser contains no transform code. Every value comes from an instrumented build of OpenFHE 1.5.1, so what you watch is what the library does.
*   **Small parameters, real code.** At $N = 8$ and $q = 17$, you can check every multiplication with a pen. The source lines are the same lines that run at production ring dimensions.

***

## References & Further Reading

If you want the mathematics and the proofs behind the code, read these next to the simulator:

1.  **The NTT algorithms that OpenFHE implements:** <a id="ref-ln"></a>
    P. Longa and M. Naehrig, *"Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography"*. CANS 2016.
    *   [Read on IACR ePrint (2016/504)](https://eprint.iacr.org/2016/504)
    *   OpenFHE cites this paper in [`transformnat.h`](https://github.com/openfheorg/openfhe-development/blob/v1.5.1/src/core/include/math/hal/intnat/transformnat.h). The code that the simulator replays implements Algorithms 1 and 2 of this paper.

2.  **The Shoup modular multiplication:** <a id="ref-harvey"></a>
    D. Harvey, *"Faster arithmetic for number-theoretic transforms"*. Journal of Symbolic Computation, 2014.
    *   [Read on arXiv (1205.2926)](https://arxiv.org/abs/1205.2926)

3.  **The scaling optimization in the inverse transform:**
    [OpenFHE issue 872](https://github.com/openfheorg/openfhe-development/issues/872) records why the division by $N$ hides inside a twiddle.

**Links and versions:**

*   The simulator: [NTT / INTT Explorer](https://visuallearn.github.io/fhe-ntt-sim/)
*   The source: [github.com/visuallearn/fhe-ntt-sim](https://github.com/visuallearn/fhe-ntt-sim)
*   OpenFHE 1.5.1, tag [`v1.5.1`](https://github.com/openfheorg/openfhe-development/tree/v1.5.1), commit [`1306d14f8c26`](https://github.com/openfheorg/openfhe-development/commit/1306d14f8c26)
*   The algorithms: [`src/core/include/math/hal/intnat/transformnat-impl.h`](https://github.com/openfheorg/openfhe-development/blob/v1.5.1/src/core/include/math/hal/intnat/transformnat-impl.h), [lines 302 to 374](https://github.com/openfheorg/openfhe-development/blob/v1.5.1/src/core/include/math/hal/intnat/transformnat-impl.h#L302-L374) for the forward transform and [lines 511 to 625](https://github.com/openfheorg/openfhe-development/blob/v1.5.1/src/core/include/math/hal/intnat/transformnat-impl.h#L511-L625) for the inverse

***

## Suggested Citation

If you found this article useful and wish to cite it in your work, we suggest:

```text
Ahmad Al Badawi, A Step-by-Step Simulator for the NTT in FHE, 2026, https://ahmadalbadawi.com/posts/2026/08/ntt-simulator-fhe/
```

Or in BibTeX:

```bibtex
@misc{nttsim2026,
  author = {Al Badawi, Ahmad},
  title = {A Step-by-Step Simulator for the NTT in FHE},
  year = {2026},
  month = {August},
  howpublished = {\url{https://ahmadalbadawi.com/posts/2026/08/ntt-simulator-fhe/}},
  note = {Accessed: [Current Date]}
}
```

***

> **Feedback:**
> Please direct any typos, questions, comments, or issues to me at my [contact page](/contact/).

***

> **License:**
> The code in this article is licensed under the **[MIT License](https://opensource.org/licenses/MIT)**.
> The text and content are licensed under **[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)**.

***
