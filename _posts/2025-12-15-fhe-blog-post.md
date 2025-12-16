---
title: 'SIMD Packing in BGV/BFV FHE Schemes'
date: 2025-12-15
permalink: /posts/2025/12/simd-packing-bgv-bfv-fhe/
tags:
  - FHE
  - BGV
  - BFV
  - SIMD Packing
  - Data Packing
  - Integer Encoding
  - Privacy
  - Homomorphic Encryption
---

***

This article provides a deep dive into the SIMD (Single Instruction, Multiple Data) packing techniques used in BGV and BFV Fully Homomorphic Encryption (FHE) schemes. We explore the mathematical isomorphism between vectors and polynomials that enables parallel processing of encrypted data. The roadmap for this article includes a high-level intuitive explanation using standard algebra, a detailed breakdown of the Chinese Remainder Theorem (CRT) as applied to cyclotomic rings, and finally, a concrete step-by-step code walkthrough in Python demonstrating data encoding and decoding procedures in BGV/BFV.

***Disclaimer**: This article is written for readers with basic background in cryptography or advanced mathematics and FHE. It assumes a working knowledge of abstract algebra concepts; specifically rings, finite fields, and polynomial arithmetic, as well as familiarity with the basic mechanics of homomorphic encryption.*

***

**Contents**
* ToC
{:toc}


***

### Introduction

Broadly speaking, FHE schemes fall into two main families: those optimized for boolean operations or low-bit-width integers, and those designed for computations on long vectors of data (e.g., real numbers, complex numbers, or integers). In this article, we focus on the BGV and BFV FHE schemes, which are particularly well-suited for processing long vectors of integers.

**SIMD parallelism is the cornerstone of these arithmetic schemes.** While modern FHE libraries abstract away the complexity of ciphertext packing, a true grasp of the technology requires looking beyond the interface. This article is for the curious FHE practitioner who wants to peel back that abstraction to understand the elegant know-how behind the machinery.

If you have used these schemes before, you know that their "magic" lies in their ability to process long vectors of integers simultaneously. A single addition or multiplication performed on a ciphertext translates to the same operation being applied to every slot of the underlying vector. This capability is what enables high-performance FHE; without it, we would be limited to encrypting integers one by one, which is prohibitively slow for complex applications.

***

### Vectors vs. Polynomials

To understand what is happening under the hood, let's look at the arithmetic interface these schemes provide. If you have used BGV ([BGV et al., 2011](#ref-bgv)) or BFV ([BFV et al., 2012](#ref-bfv)), you are essentially given machinery to perform operations on encrypted integer vectors.

Formally, given input vectors $\mathbf{u}$ and $\mathbf{v}$, the encryption process creates objects that satisfy specific homomorphic properties:

$$
\text{ctxt}_{\mathbf{u}} = \text{Encrypt}(\text{Encode}(\{u_0, u_1, \ldots, u_{n-1}\}))
$$

$$
\text{ctxt}_{\mathbf{v}} = \text{Encrypt}(\text{Encode}(\{v_0, v_1, \ldots, v_{n-1}\}))
$$

The schemes are designed such that adding these ciphertexts yields the encryption of the component-wise sum:

$$
\text{ctxt}_{\text{sum}} = \text{ctxt}_{\mathbf{u}} + \text{ctxt}_{\mathbf{v}} = \text{Encrypt}(\text{Encode}(\{u_0 + v_0, \ldots, u_{n-1} + v_{n-1}\}))
$$

Similarly, multiplication yields the encryption of the component-wise product:

$$
\text{ctxt}_{\text{prod}} = \text{ctxt}_{\mathbf{u}} \times \text{ctxt}_{\mathbf{v}} = \text{Encrypt}(\text{Encode}(\{u_0 \times v_0, \ldots, u_{n-1} \times v_{n-1}\}))
$$

This property, where one ciphertext operation equals $n$ plaintext operations, is the **SIMD** breakthrough pioneered by Smart and Vercauteren ([Smart & Vercauteren, 2011](#ref-simd)).

**However, a structural disconnect exists.** Internally, these FHE schemes do not operate on vectors; they operate on large polynomials. This raises a critical question: **How can a polynomial efficiently represent and process a vector of independent data?** The answer lies in the **Chinese Remainder Theorem (CRT)**, which provides the mathematical isomorphism necessary to pack a vector of numbers into a single polynomial.

This article will unpack the core concept of SIMD packing, the engine behind popular ring-based HE schemes, using a concrete algebraic example. You can follow along with the code used in this article in the **[accompanying GitHub repository](https://github.com/caesaretos/bgv_bfv_slots_packing)**.

***

## Problem: FHE Requires Polynomials, But We Need Vectors

To understand how FHE enables vector operations, it is crucial to distinguish between **encoding** and **encryption**. While encryption provides the security (the "lock on the box"), encoding provides the mathematical structure that makes computation possible (the "internal machinery"). The 'magic' of SIMD parallelism is dictated by the algebraic properties of the encoding layer, the Plaintext Space, not the encryption itself. The encryption process merely adds cryptographic noise to this structure without altering its fundamental arithmetic properties.

For this reason, we will temporarily disregard the encryption step. If the desired arithmetic operation holds true for the encoded plaintexts, the defining homomorphic property of FHE guarantees that the same operation will yield the correct result when performed on the encrypted data as well. This is the mathematical essence of a homomorphism: operations conducted in one space are faithfully preserved when translated into a structurally similar second space. Our focus is exclusively on the encoding map, as this is where vector operations are mathematically defined.

Considering only the encoding layer, we can simplify Equations 1-4 above as follows:

$$\text{ptxt}_{\mathbf{u}} = \text{Encode}(\{u_0, u_1, \ldots, u_{n-1}\})$$  

$$\text{ptxt}_{\mathbf{v}} = \text{Encode}(\{v_0, v_1, \ldots, v_{n-1}\})$$ 

Adding these plaintexts results in point-wise addition of the operand vectors as follows:

$$\text{ptxt}_{\text{sum}} = \text{ptxt}_{\mathbf{u}} + \text{ptxt}_{\mathbf{v}}$$

Decoding to expected result:

$$\text{Decode}(\text{ptxt}_{\text{sum}}) = \{u_0 + v_0, u_1 + v_1, \ldots, u_{n-1} + v_{n-1}\}$$

Similarly, multiplying these plaintexts results in point-wise multiplication of the operand vectors:

$$\text{ptxt}_{\text{prod}} = \text{ptxt}_{\mathbf{u}} \times \text{ptxt}_{\mathbf{v}}$$

Decoding to expected result:

$$\text{Decode}(\text{ptxt}_{\text{prod}}) = \{u_0 \times v_0, u_1 \times v_1, \ldots, u_{n-1} \times v_{n-1}\}$$

In BGV/BFV, a single plaintext is represented by a polynomial within a specialized mathematical structure called a ring. This polynomial acts as a container, or 'slot' system, holding an array of multiple data values. For this reason, polynomials are the core mathematical object, and all homomorphic calculations rely exclusively on polynomial arithmetic.

We have established that the fundamental data unit in BGV and BFV schemes is the polynomial; consequently, all cryptographic operations are performed as polynomial arithmetic. This leads us to the central question of this article: 

<div style="background-color: #f8f9fa; border-left: 4px solid #333; padding: 15px; margin: 20px 0;" markdown="1">
How do you map data vectors $$\mathbf{u} = \{u_0, u_1, \ldots, u_{n-1}\}$$ and $$\mathbf{v} = \{v_0, v_1, \ldots, v_{n-1}\}$$ into polynomials $U(x)$ and $V(x)$ such that, for an operation $$\star \in \{+, \times\}$$, the resulting polynomial $U(x) \star V(x)$ satisfies the following decoding property:

$$ \text{Decode}(U(x) \star V(x)) = \{u_0 \star v_0, u_1 \star v_1, \ldots, u_{n-1} \star v_{n-1}\} $$
</div>

***
## Basic Intuition: How Polynomials Affect Vector Operations in a Single Shot

Before we dive into complex cryptography or FHE, we need to appreciate a fundamental mathematical property of polynomials.

In standard programming, if you want to add two arrays of numbers, you use a loop to add them index by index. But in algebra, we can use **polynomial interpolation** to treat an entire array of numbers as a single mathematical object. By performing one operation on the polynomial, we effectively perform that operation on every number inside it simultaneously.

Here is a simple demonstration of how this works, using standard arithmetic. Although this example accurately illustrates the core intuition behind SIMD encoding, we note that the specific underlying algebraic rings used in BGV/BFV schemes are more complex than those demonstrated in this simple example. The key takeaway, however, is the concept of leveraging polynomial arithmetic for vector operations.

### 1. Setup: Data as Coordinates
Imagine we have two data vectors, $\mathbf{u}$ and $\mathbf{v}$. To transform them into polynomials, we set up a coordinate system by pre-selecting a fixed list of $x$-coordinates. The data values from the vectors are then treated as the $y$-coordinates, creating a set of fixed points for interpolation.

*   **Vector $\mathbf{u}$:** $$\{1, 10, 49, 142\}$$ at positions $$\mathbf{x} = \{0, 1, 2, 3\}$$
*   **Vector $\mathbf{v}$:** $$\{5, 2, 5, 32\}$$ at positions $$\mathbf{x} = \{0, 1, 2, 3\}$$  
Note that the positions ($\mathbf{x}$) must be the same for each vector we want to process, but we are free to choose the positions.

We can define these data points in Python using the `sympy` library to handle the symbolic math:

```python
import sympy as sp
x = sp.Symbol('x')

# 1. Define Data Points
data_u = [(0, 1), (1, 10), (2, 49), (3, 142)]
data_v = [(0, 5), (1, 2), (2, 5), (3, 32)]
```

### 2. Data Encoding: Interpolation
Instead of storing the vectors as lists, we find a curve (a polynomial) that passes perfectly through these points.
*   **$U(x)$** is the equation that passes through every point in Vector $\mathbf{u}$.
*   **$V(x)$** is the equation that passes through every point in Vector $\mathbf{v}$.

In our example, $U(x) = 1+2x+3x^2+4x^3$ and $V(x) = 5-6x^2+3x^3$. You can verify that by evaluating $U(x)$ (resp. $V(x)$) at the positions and that should give you $\mathbf{u}$ (resp. $\mathbf{v}$). The required interpolating polynomial can be found efficiently using techniques like Lagrange interpolation.

Here is how we perform this interpolation programmatically to obtain our algebraic objects:

```python
# 2. Interpolate the Polynomials over Q (the set of rational numbers)
U_x = sp.interpolate(data_u, x)
V_x = sp.interpolate(data_v, x)

print("\n Interpolated Polynomials")
print("U(x):", sp.expand(U_x))
print("V(x):", sp.expand(V_x))
```
**Output:**
```text
Interpolated Polynomials
U(x): 4*x**3 + 3*x**2 + 2*x + 1
V(x): 3*x**3 - 6*x**2 + 5
```

### 3. "Single Shot" Operations
Now comes the magic. We don't touch the numbers inside the vectors. Instead, we simply add or multiply the **polynomial equations themselves**.

$$ \text{Sum}(x) = U(x) + V(x) = 7x^3 - 3x^2 + 2x + 6 $$

$$ \text{Prod}(x) = U(x) \times V(x) = 12x^6 - 15x^5 - 12x^4 + 11x^3 + 9x^2 + 10x + 5 $$

In our Python simulation, because we are working with symbolic objects, performing these component-wise operations on the underlying data is as simple as adding or multiplying the polynomial variables:

```python
# 3. Compute Sum and Product Polynomials (using schoolbook algebra)
Sum_x = sp.expand(U_x + V_x)
Prod_x = sp.expand(U_x * V_x)

print("\n--- Resulting Polynomials ---")
print("Sum(x):", Sum_x)
print("Prod(x):", Prod_x)
```
**Output:**
```text
--- Resulting Polynomials ---
Sum(x): 7*x**3 - 3*x**2 + 2*x + 6
Prod(x): 12*x**6 - 15*x**5 - 12*x**4 + 11*x**3 + 9*x**2 + 10*x + 5
```

### 4. Verification
Does this algebraic manipulation actually affect the data points correctly? Let's check the result by evaluating our new polynomials at the original positions ($$\mathbf{x}=\{0, 1, 2, 3\}$$).

| Index ($\mathbf{x}$) | $\mathbf{u}$ | $\mathbf{v}$ | **$\text{Sum}(x)$** | **Expected Sum** | **$\text{Prod}(x)$** | **Expected Product** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **0** | 1 | 5 | **6** | $1+5=6$ | **5** | $1\times5=5$ |
| **1** | 10 | 2 | **12** | $10+2=12$ | **20** | $10\times2=20$ |
| **2** | 49 | 5 | **54** | $49+5=54$ | **245** | $49\times5=245$ |
| **3** | 142 | 32 | **174** | $142+32=174$ | **4544** | $142\times32=4544$ |

We can run a verification loop in our script to demonstrate that evaluating the resulting polynomials matches the manual arithmetic:

```python
# 4. Verification
print("\n--- Verification of Component-Wise Operation ---")
print("X-Value | U(x) | V(x) | Sum(x) | U(x)+V(x) | Prod(x) | U(x)*V(x)")
print("-" * 100)

for i in [0, 1, 2, 3]:
    p1_val = U_x.subs(x, i)
    p2_val = V_x.subs(x, i)
    psum_val = Sum_x.subs(x, i)
    pprod_val = Prod_x.subs(x, i)

    # Convert everything to string first for formatting
    i_s = str(i)
    p1_s = str(p1_val)
    p2_s = str(p2_val)
    psum_s = str(psum_val)
    ysum_s = str(p1_val + p2_val)
    pprod_s = str(pprod_val)
    yprod_s = str(p1_val * p2_val)

    print(f"{i_s:^7} | {p1_s:^4} | {p2_s:^4} | {psum_s:^6} | {ysum_s:^9} | {pprod_s:^7} | {yprod_s:^6}")
```
**Output:**
```text
--- Verification of Component-Wise Operation ---
X-Value | U(x) | V(x) | Sum(x) | U(x)+V(x) | Prod(x) | U(x)*V(x)
----------------------------------------------------------------------------------------------------
   0    |  1   |  5   |   6    |     6     |    5    |   5   
   1    |  10  |  2   |   12   |    12     |   20    |   20  
   2    |  49  |  5   |   54   |    54     |   245   |  245  
   3    | 142  |  32  |  174   |    174    |  4544   |  4544 
```

This demonstration confirms the validity of the SIMD arithmetic approach. Specifically, the evaluation of the resulting sum or product polynomial at each pre-selected position (slot) precisely matches the component-wise arithmetic operation applied to the original input vectors. This demonstrates that algebraic manipulations of a single encoded polynomial can simultaneously process all the underlying data slots.

### 5. Key Takeaway
As you can see, adding the polynomials resulted in a new curve that represents the **sum** of the data at every position point. Multiplying them created a curve representing the **product**.

We successfully processed four independent data points using a single algebraic instruction on polynomials. This ability to "pack" data into a polynomial and process it in parallel is the mathematical engine that powers modern FHE schemes.

***

## Back to FHE: Chinese Remainder Theorem and Polynomial Roots

While the example above effectively illustrated the intuition behind polynomial packing, it glossed over critical engineering constraints. Real-world FHE schemes like BGV and BFV cannot operate over simple rational numbers ($\mathbb{Q}$) because they lack bounds. You may have noticed that the coefficients in our example grew rapidly; in a finite computer system, this unbounded growth leads to precision overflows. Furthermore, standard polynomial multiplication causes the degree to increase (e.g., multiplying two degree-$3$ polynomials results in a degree-$6$ polynomial). To prevent data from growing unboundedly in both size and dimension, BGV and BFV schemes employ **Modular Arithmetic** and **Cyclotomic Quotient Rings**. These structures impose strict mathematical limits, ensuring that both the coefficient size and the polynomial degree remain fixed, regardless of how many operations we perform.

In these schemes, a single plaintext is a polynomial $U(x)$ residing in a specific ring, typically $\mathcal{R} = \mathbb{Z}_p[x]/\langle \Phi_m(x) \rangle$, where $\Phi_m(x)$ is a cyclotomic polynomial (often $x^{n}+1$). We assume here that both $m$ and $n$ are powers of two.

Think of the objects within this ring as regular polynomials but with degree less than $n$ and their coefficients are integers modulo $p$, that is, $\in \mathbb{Z}_p$. This raises the question: **What serves as the "coordinate system" for our interpolation in this modular world?**

The answer lies in the **Chinese Remainder Theorem (CRT)**.

### 1. The Ring and the Roots (The "Slots")

In these FHE schemes, the mathematical ring is carefully selected to support SIMD operations. We specifically choose the modulus $p$ such that the polynomial $x^n+1$ splits into $n$ distinct linear factors. This occurs when $p \equiv 1 \pmod{m}$ (where $m=2n$).

Under this condition, there exist $n$ distinct roots (values) $$\{r_0, r_1, \ldots, r_{n-1}\}$$ in $\mathbb{Z}_p$ such that:

$$
x^n + 1 \equiv (x - r_0)(x - r_1)\cdots(x - r_{n-1}) \pmod{p}
$$

In our previous Python example, we arbitrarily chose the positions $$\{0, 1, 2, 3\}$$ as our $x$-coordinates. In BGV/BFV schemes, these roots $r_i$ become our fixed $x$-coordinates, serving as the **encoding slots**. In BGV/BFV jargon, these roots are formally known as primitive $m$-th roots of unity (for some power-of-two integer $m$). There are always exactly $n$ such roots that can be used as the discrete positions for the slots.

### 2. The CRT Bridge: Decomposition vs. Composition

To intuitively understand why we use the **Inverse** CRT for packing, it helps to view the relationship between polynomials and vectors through the lens of **Decomposition** and **Synthesis**.

*   **The Polynomial** ($U(x)$) is the **"Macroscopic Entity."** It is a single, monolithic object that the FHE scheme operates on directly.
*   **The Slots** (the values at roots $r_i$) are the **"Microscopic Entities."** These are the independent, smaller components that constitute the whole.

The standard **Forward CRT** acts as a mechanism for **decomposition**: it takes the big entity (the polynomial) and breaks it down into its constituent parts (the values at the roots).

$$ \text{Forward CRT (Decomposition): } \text{Polynomial } \rightarrow \text{Vector of Parts} $$

However, our scenario presents the reverse problem. We start with the vector of data, we already possess the discrete, independent entities, and we need to combine them into a single, cohesive structure that the FHE scheme can encrypt and operate on. We are not analyzing a polynomial; we are building one.

Therefore, we must travel in the opposite direction: from the parts to the whole. This process of **Synthesis** is exactly what the **Inverse CRT** (Interpolation) achieves.

$$ \text{Inverse CRT (Synthesis): } \text{Vector of Parts } \rightarrow \text{Polynomial} $$

***

### 3. Packing (Encoding) via Inverse CRT (Interpolation)

The CRT guarantees that a polynomial $U(x) \in \mathcal{R}$ is uniquely determined by its values when evaluated at these $n$ distinct roots. This gives us a direct mapping between our data vector and the polynomial:

$$
\text{Vector } \mathbf{u} = \{u_0, \ldots, u_{n-1}\} \leftrightarrow \text{Polynomial } U(x)
$$

This mapping relies on the condition that $U(x=r_i) \equiv u_i \pmod{p}$ for all $i=0$ to $n-1$.

To pack the vector $\mathbf{u}$, we perform the same **interpolation** process (Inverse CRT) we saw earlier, but this time using the pairs of (root, data value): $\langle (r_0, u_0), (r_1, u_1), \ldots, (r_{n-1}, u_{n-1}) \rangle$. 

A natural question arises regarding the specific sequence of the roots $r_i$. While the ordering of roots is critical for implementing efficient **rotation operations** (cyclic shifts) on encrypted vectors, it does not affect the element-wise addition or multiplication discussed here. Therefore, we will assume a fixed arbitrary order for this article and reserve the intricacies of root permutation and Galois automorphisms for a future article.

### 4. SIMD Operation on Ciphertext

This is where the algebraic structure pays off. Operations on the polynomials in the ring $\mathcal{R}$ translate *directly* into component-wise operations on the vectors encoded within them.

If $U(x)$ encodes vector $\mathbf{u}$ and $V(x)$ encodes vector $\mathbf{v}$, then:

$$
U(x) + V(x) \text{ encodes } \{u_0+v_0, u_1+v_1, \ldots, u_{n-1}+v_{n-1}\}
$$

and

$$
U(x) \times V(x) \text{ encodes } \{u_0\times v_0, u_1\times v_1, \ldots, u_{n-1}\times v_{n-1}\}
$$

This works because for any root $r_i$, the evaluation is a homomorphism:

$$
(U(x) + V(x)) \big|_{x=r_i} = U(r_i) + V(r_i) \equiv u_i + v_i \pmod{p}
$$

and

$$
(U(x) \times V(x)) \big|_{x=r_i} = U(r_i) \times V(r_i) \equiv u_i \times v_i \pmod{p}
$$

**A single polynomial addition achieves $n$ independent vector additions simultaneously.**

### 5. Unpacking (Decoding) via Forward CRT

To recover the resulting data vector from a result polynomial $\text{Sum}(x) = U(x) + V(x)$ or $\text{Prod}(x) = U(x) \times V(x)$, we simply evaluate $\text{Sum}(x)$ or $\text{Prod}(x)$ at all $n$ roots:

$$
\text{Decoded Vector} = \{\text{Sum}(r_0), \text{Sum}(r_1), \ldots, \text{Sum}(r_{n-1})\}
$$

or

$$
\text{Decoded Vector} = \{\text{Prod}(r_0), \text{Prod}(r_1), \ldots, \text{Prod}(r_{n-1})\}
$$

***

## Code Walkthrough: A Concrete Algebraic Example

To visualize this with actual modular arithmetic, let's look at a concrete setup using Python and the **SymPy** library, which allows us to manipulate polynomials over finite fields easily.

### Setup

| Parameter | Value | Description |
| :--- | :--- | :--- |
| **$m$** | $8$ | Cyclotomic Index (defines the field) |
| **$p$** | $73$ | Prime Modulus ($73 \equiv 1 \pmod{8}$) |
| **$n$** | $4$ | Polynomial Degree ($n = m/2$) |
| **Ring** | $\mathbb{Z}_{73}[x] / \langle x^4 + 1 \rangle$ | The mathematical space |

First, we set up our environment and find the $n=4$ distinct roots (the slots) in $\mathbb{Z}_{73}$ where $x^4+1 \equiv 0 \pmod{73}$.

```python
import sympy
from sympy import symbols, Poly, GF

# 1. SETUP THE RING PARAMETERS
m = 8
p = 73
n = m // 2

x = symbols('x')
# Define the ring modulus polynomial: x^n + 1
ring_modulus = Poly(x**n + 1, x, domain=GF(p))

# Find the "Slots" (roots of unity) in a very lazy way!
roots = [r for r in range(p) if (pow(r, n, p) + 1) % p == 0]

print(f"Setup: Ring Z_{p}[x] / (x^{n} + 1)")
print(f"Found {len(roots)} slots (roots) at indices: {roots}")
```
**Output:**
```text
--------------------------------------------------
Setup: Ring Z_73[x] / (x^4 + 1)
Found 4 slots (roots) at indices: [10, 22, 51, 63]
```

### Encoding (Packing)

We define two data vectors, $\mathbf{u}$ and $\mathbf{v}$, that we wish to process.

$$ \mathbf{u} = \{0, 1, 2, 3\} $$ and $$ \mathbf{v} = \{4, 5, 6, 7\} $$

To pack these into polynomials, we perform Lagrange Interpolation over the finite field. We use a helper function, `inverse_crt_interpolate`, which constructs the unique polynomial that evaluates to our data values at the specified roots. *(Note: The full implementation of this interpolation function is available in the associated source code on GitHub)*.

```python
# 2. DEFINE INPUT DATA
# Ensure inputs are modulo p
u = [0, 1, 2, 3]
v = [4, 5, 6, 7]

# ==========================================
# 3. PACKING (ENCODING)
# ==========================================
points_u = list(zip(roots, u))
points_v = list(zip(roots, v))

# Perform Interpolation using our custom Z_p function
packed_poly_U = inverse_crt_interpolate(points_u, x, p)
packed_poly_V = inverse_crt_interpolate(points_v, x, p)

print("-" * 50)
print("PACKING STEP")
print("Polynomial U(x):")
print(packed_poly_U.as_expr())
print("\nPolynomial V(x):")
print(packed_poly_V.as_expr())
```
**Output:**
```text
--------------------------------------------------
PACKING STEP
Polynomial U(x):
13*x**3 - 19*x - 35

Polynomial V(x):
13*x**3 - 19*x - 31
```

### SIMD Check

Now we perform the arithmetic operations directly on the polynomials.
1.  **Addition:** We simply add the polynomials.
2.  **Multiplication:** We simply multiply the polynomials.

Both additions and multiplications must be carried out in the plaintext ring, that is, we take the remainder modulo $$\{(x^4+1), p\}$$ to ensure the result stays within the ring. This modular reduction is what ensures our resultant polynomials are bounded in both coefficient size and degree.

```python
# ==========================================
# 4. SIMD OPERATION CHECK
# ==========================================
# Addition: Element-wise addition in the ring
sum_poly = (packed_poly_U + packed_poly_V) % ring_modulus

# Multiplication: Element-wise multiplication in the ring
# IMPORTANT: In the ring Z_p[x] / (x^n + 1), multiplication increases degree.
# We must reduce the result modulo (x^n + 1).
prod_poly = (packed_poly_U * packed_poly_V) % ring_modulus

print("-" * 50)
print("Polynomial Sum(x):")
print(sum_poly.as_expr())
print("\nPolynomial Prod(x):")
print(prod_poly.as_expr())
```
**Output:**
```text
--------------------------------------------------
Polynomial Sum(x):
26*x**3 + 35*x + 7

Polynomial Prod(x):
18*x**3 - 27*x**2 + 13*x - 27
```

### Decoding (Unpacking)

To verify the result, we evaluate the resulting `sum_poly` and `prod_poly` at the original roots. If the SIMD property holds, these evaluations should match the component-wise sum and product of our original vectors $\mathbf{u}$ and $\mathbf{v}$.

```python
# ==========================================
# 5. UNPACKING (DECODING)
# ==========================================
# Evaluate at roots.
# .eval() returns a GF element, we cast to int for display.
recovered_data_U = [int(packed_poly_U.eval(r)) for r in roots]
recovered_data_V = [int(packed_poly_V.eval(r)) for r in roots]
sum_data = [int(sum_poly.eval(r)) for r in roots]
prod_data = [int(prod_poly.eval(r)) for r in roots]

print("-" * 50)
print("UNPACKING STEP")

# Check Inputs
is_match_u = (recovered_data_U == u)
print(f"Recovered U matches input? {is_match_u}")

is_match_v = (recovered_data_V == v)
print(f"Recovered V matches input? {is_match_v}")

# Check Sum
print("-" * 50)
print("SIMD CHECK (U + V)")
target_sum = [(i + j) % p for i, j in zip(u, v)]
print(f"Expected Sum: {target_sum}")
print(f"Actual Sum:   {sum_data}")
print(f"Match?        {sum_data == target_sum}")

# Check Product
print("-" * 50)
print("SIMD CHECK (U * V)")
target_prod = [(i * j) % p for i, j in zip(u, v)]
print(f"Expected Prod: {target_prod}")
print(f"Actual Prod:   {prod_data}")
print(f"Match?         {prod_data == target_prod}")
```
**Output:**
```text
--------------------------------------------------
UNPACKING STEP
Recovered U matches input? True
Recovered V matches input? True
--------------------------------------------------
SIMD CHECK (U + V)
Expected Sum: [4, 6, 8, 10]
Actual Sum:   [4, 6, 8, 10]
Match?        True
--------------------------------------------------
SIMD CHECK (U * V)
Expected Prod: [0, 5, 12, 21]
Actual Prod:   [0, 5, 12, 21]
Match?         True
```

As expected, the operations on the polynomials resulted in every element of the packed vectors being successfully added and multiplied simultaneously!

***

## Key Takeaways

*   **SIMD is Essential for Throughput:**
    The ability to perform **SIMD** operations via polynomial encoding is what makes Homomorphic Encryption schemes computationally practical. Without SIMD, you would need one massive polynomial ciphertext for every single number you wanted to process. With SIMD, you process thousands of data points for the "price" of one.

*   **The Modulus $p$ Defines Data Type/Bit-Width:**
    Just as a standard `int8` can only hold values up to 127, the prime modulus $p$ acts as the hard ceiling for your data. It dictates the **bit-width** and precision of the integers you can compute on. If your calculation (e.g., summing many vectors) exceeds $p$, the values will "wrap around" modulo $p$. If your application expects integer arithmetic over $\mathbb{Z}$, this wraparound will yield incorrect results. Therefore, $p$ must be chosen large enough to accommodate the maximum possible result of your computation.

*   **The Degree $n$ Defines Batch Size:**
    The polynomial degree $n$ (derived from the cyclotomic index $m$) determines the **capacity** of your ciphertext. Specifically, you get exactly $n$ independent slots. Larger $n$ means more parallelism (higher throughput), but it also results in larger ciphertexts and slower individual polynomial operations.

*   **CRT is the "Switching Engine":**
    The Chinese Remainder Theorem isn't just a math trick; it is the translation engine. It provides the **isomorphism** (a 1-to-1 structural mapping) that allows us to switch between the "Polynomial Domain" (where FHE works efficiently) and the "Vector Domain" (where the data makes sense to us).

*   **Amortized Efficiency:**
    Because the FHE scheme operates on the polynomial as a whole, the computational cost of adding two ciphertexts is the same regardless of whether they contain 1 packed item or 10,000 packed items, assuming a fixed polynomial degree. This makes FHE uniquely suited for "wide" workloads like matrix multiplication or database queries.

This technique is a perfect example of how abstract algebra can be translated into a powerful tool for modern, privacy-preserving computation.

## References & Further Reading

If you want to dive deeper into the mathematics and proofs behind these concepts, here are the foundational papers that established the standard for SIMD and Ring-based BGV and BFV schemes:

1.  **The SIMD Breakthrough (Packing):** <a id="ref-simd"></a>
    N. P. Smart and F. Vercauteren, *"Fully Homomorphic SIMD Operations"*. Designs, Codes and Cryptography, 2014.
    *   [Read on IACR ePrint (2011/133)](https://eprint.iacr.org/2011/133)

2.  **The BGV Scheme:** <a id="ref-bgv"></a>
    Z. Brakerski, C. Gentry, and V. Vaikuntanathan, *"Fully Homomorphic Encryption without Bootstrapping"*. TOCT 2014.
    *   [Read on IACR ePrint (2011/277)](https://eprint.iacr.org/2011/277)

3.  **The BFV Scheme:** <a id="ref-bfv"></a>
    J. Fan and F. Vercauteren, *"Somewhat Practical Fully Homomorphic Encryption"*. IACR ePrint 2012.
    *   [Read on IACR ePrint (2012/144)](https://eprint.iacr.org/2012/144)

***

## Citing this Article

If you found this article useful and wish to cite it in your work, please use the following BibTeX entry:

```bibtex
@misc{simdpacking2025,
  author = {Al Badawi, Ahmad},
  title = {SIMD Packing in BGV/BFV FHE Schemes},
  year = {2025},
  month = {December},
  howpublished = {\url{https://ahmadalbadawi.com/posts/2025/12/simd-packing-bgv-bfv-fhe/}},
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
