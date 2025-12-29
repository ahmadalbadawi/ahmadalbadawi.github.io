---
title: 'The Math of Slot Rotation in BGV and BFV'
date: 2025-12-28
permalink: /posts/2025/12/slot-rotation-bgv-bfv-fhe/
tags:
  - FHE
  - BGV
  - BFV
  - SIMD Packing
  - Slot Rotation
  - Automorphism
  - Integer Encoding
  - Privacy
  - Homomorphic Encryption
---

**[Ahmad Al Badawi](https://ahmadalbadawi.com/)**

***

This article extends the previous article on [SIMD Packing in BGV/BFV FHE Schemes](https://ahmadalbadawi.com/posts/2025/12/simd-packing-bgv-bfv-fhe/), where we explored how to encode a vector of messages into a single polynomial using the roots of unity. While that article demonstrated how to encode and decode data, and to do basic arithmetic such as point-wise addition and multiplication efficiently, it left us with a static object (a polynomial): the values were rooted in their respective slots without a method for structured shuffling or permutation.

While point-wise addition and multiplication are powerful, they are insufficient for performing arbitrary computations on encoded vectors. Because these operations are strictly component-wise, data at different indices remain isolated. In the previous article, we noted that this functionality did not depend on the ordering of the data slots. To break this isolation and enable arbitrary logic, we must now address this very topic: by choosing a *specific*, non-arbitrary ordering of the slots, we can implement permutation operators (such as cyclic shifts) to facilitate interaction across the entire vector.

In this article, we will explore how vector permutations can be achieved via simple polynomial operations. We will dive into **Galois Theory** to understand how the abstract "symmetry" of roots allows us to shuffle data, and how we can utilize **Cyclic Groups** to turn those shuffles into precise, controllable rotations. It was the efficient use of these permutations that allowed Gentry, Halevi, and Smart ([GHS et al., 2012](#ref-simd-auto)) to demonstrate that FHE has only "polylogarithmic overhead", proving that encrypted computation does not need to be exponentially slower than plaintext computation.

Just like our previous article, we will disregard the encryption layer to focus on the encoding map. This is where the vector operations are actually defined. Thanks to the homomorphic property of FHE schemes, we know that whatever logic we establish for encoding will automatically hold true for the encrypted ciphertext. The encryption layer simply maintains data confidentiality across all performed operations. 

You can follow along with the code used in this article in the **[accompanying GitHub repository](https://github.com/caesaretos/bgv_bfv_slots_packing)**.

***Disclaimer**: This article is written for readers with basic background in cryptography or advanced mathematics and FHE. It assumes a working knowledge of abstract algebra concepts; specifically rings, finite fields, and polynomial arithmetic, as well as familiarity with the basic mechanics of homomorphic encryption.*

***

**Contents**
* ToC
{:toc}


***

## Introduction

In our previous discussion on [SIMD packing in BFV and BGV](https://ahmadalbadawi.com/posts/2025/12/simd-packing-bgv-bfv-fhe/), we learned how to take a vector of individual values, say $\mathbf{v} = (v_0, v_1, \dots, v_{n-1})$, and encode it into a single polynomial. This was achieved by constructing a special polynomial $V(x)$ such that evaluating $V(x)$ at specific, predetermined roots of unity ($r_0, r_1, \dots$) would reveal our individual messages ($V(x=r_i) = v_i$). This technique is incredibly powerful, allowing us to encode and operate on multiple values simultaneously, significantly boosting efficiency in FHE schemes like BGV ([BGV et al., 2011](#ref-bgv)) and BFV ([BFV et al., 2012](#ref-bfv)). We focused on component-wise addition and multiplication, which are achieved by performing those same operations on the polynomials encoding the original vectors.

However, simply packing data is not enough for many real-world applications. Imagine performing a convolution (a fundamental operation in signal processing and deep learning) or a matrix multiplication. These tasks inherently require shifting elements of a vector, e.g., moving $v_0$ to the position where $v_1$ used to be, and so on. If our data remains static and rooted in its original slots within the polynomial, these advanced computations become very challenging. The main challenge is to perform these shifts, or "rotations", on the polynomials encoding the original vectors using polynomial operations only.

So, how do we achieve this seemingly impossible feat? The answer lies not in directly manipulating the seemingly random coefficients of the encoding polynomial but in a deeper mathematical concept. We do not try to "shift" the polynomial directly. Instead, we leverage the elegant properties of **Galois Automorphisms**. These powerful operators allow us to systematically permute the underlying roots of unity, and by doing so, subtly rearrange the messages stored in their associated slots, all while maintaining the integrity of data encoding. This transformation is the key to provide arbitrary computations in the FHE world.

***

## Galois Theory and Symmetry

Before we dive deep into how to rotate a vector through its encoding polynomial, we need to understand *why* we are allowed to manipulate the encoding polynomial in this specific way. The answer lies in a branch of mathematics called **Galois Theory**, which studies the "symmetry" of roots.

### The Concept: Indistinguishability (The Shell Game)
Imagine a magician performing a shell game with identical marbles. If the magician swaps two marbles while your back is turned, the game looks exactly the same to you. The marbles are indistinguishable; they have the same weight, size, and color.

In our encoding scheme, the **roots of unity** (the locations where we store our data) are like these marbles.
We are working with the polynomial $x^N + 1$ (where $N = m/2$ when $m$ is a power of two) over a finite field. The roots of this polynomial are the "Primitive $m$-th roots of unity."

To the underlying number system (modular arithmetic), these roots are **algebraically indistinguishable**.
*   They all satisfy $r_i^m \equiv 1$.
*   They all satisfy $r_i^{N} \equiv -1$.
*   There is no algebraic test we can run *inside the system* to single out one root as "primary" or "more primitive" compared to another.

### The Theory: Symmetry and the Galois Group
Because these roots are indistinguishable, we can "swap" them, and the mathematical equations that govern the encoding remain valid.

*   **The Swap:** A valid substitution of one root for another is called an **Automorphism**.
*   **The Group:** The collection of *all* possible valid swaps is called the **Galois Group**.

In the context of cyclotomic polynomials (used in BGV/BFV), the Galois Group is isomorphic to the multiplicative group of integers modulo $m$, denoted as $(\mathbb{Z}/m\mathbb{Z})^\times$. This means that every integer $k$ that is co-prime to $m$ represents a unique, valid way to shuffle the roots.

### The Mechanism: The Map $x \to x^k$
How do we actually perform this "swap" on a polynomial $U(x)$? We cannot just pick up coefficients and move them. Instead, we perform a variable substitution.

We replace every instance of $x$ with $x^k$.

$$ \sigma_k: U(x) \mapsto U(x^k) $$

Here is what happens when we do this:
1.  Our data was originally stored at a specific root, say $r$.
2.  When we transform the polynomial to $U(x^k)$, the value that *was* at $r$ moves.
3.  Effectively, the slot that used to look at root $r$ is now looks at root $r^k$.
4.  Since $r^k$ is *also* a primitive root (because $k$ is co-prime to $m$), we have not destroyed data; we have simply **permuted** which data sits in which slot.

---

### Python Illustration: Seeing the Shuffle

Let's make this concrete. We will use the parameters $m=16$ and $p=17$.
*   $m=16$: We are working with 16th roots of unity.
*   $p=17$: Our modular arithmetic is mod 17.
*   $N=8$: We have 8 slots (degree of polynomial $x^8 + 1$).

We will generate the roots, encode a simple sequence `[0, 1, ... 7]`, apply a Galois automorphism ($x \to x^3$), and see how the data moves.

*(Note: The helper functions `get_roots_of_unity`, `inverse_crt_interpolate`, and `apply_automorphism` used below are available in the **[associated GitHub repository](https://github.com/caesaretos/bgv_bfv_slots_packing)** for this article.)*

```python
from sympy import symbols

# Setup Parameters
m = 16       # Cyclotomic index
p = 17       # Modulus
N = m // 2   # Number of slots (8)
x = symbols('x')

print(f"--- SETUP: m={m}, p={p}, N={N} ---")

# A. Identify the Roots (The Slots)
# We find all primitive 16-th roots of unity modulo 17
roots = get_roots_of_unity(m, p)
roots = sorted(roots) # Sort to ensure deterministic slot ordering
print(f"Roots (Slots):        {roots}")

# B. Define Messages
msgs = [0, 1, 2, 3, 4, 5, 6, 7]
print(f"Original Message:     {msgs}")

# C. Encode (Pack Messages into Polynomial)
# We interpolate to find A(x) such that A(root_i) = msg_i
points_msgs = list(zip(roots, msgs))
poly_sym = inverse_crt_interpolate(points_msgs, x, p)

# Convert SymPy Poly to simple list of coefficients [c0, c1, c2...]
coeffs = [int(poly_sym.coeff_monomial(x**i)) % p for i in range(N)]
print(f"Encoded Poly Coeffs:  {coeffs}")
```

**Output:**
```text
--- SETUP: m=16, p=17, N=8 ---
Roots (Slots):        [3, 5, 6, 7, 10, 11, 12, 14]
Original Message:     [0, 1, 2, 3, 4, 5, 6, 7]
Encoded Poly Coeffs:  [12, 12, 0, 15, 0, 11, 0, 11]
```

Now that our data is packed, let's apply the **Galois Automorphism**. We choose $k=3$ (which is a valid action since it is co-prime to $m=16$) and map $x \to x^3$. Notice that while the **set** of roots remains identical, their **order** has been shuffled.

```python
# D. Apply Galois Automorphism (The Rotation)
k = 3
print(f"\n--- Applying Automorphism x -> x^(k={k}) ---")

# 1. Visualize how roots shuffle physically
# We calculate r^k mod p for every root r to see where they go
shuffled_roots = [pow(r, k, p) for r in roots]
print(f"Shuffled Roots:         {shuffled_roots}")

# 2. Apply automorphism to polynomial coefficients
# This computes A_new(x) = A(x^k)
rotated_coeffs = apply_automorphism(coeffs, k, N, p)
print(f"Rotated Poly Coeffs:  {rotated_coeffs}")
```

**Output:**
```text
--- Applying Automorphism x -> x^(k=3) ---
Shuffled Roots:         [10, 6, 12, 3, 14, 5, 11, 7]
Rotated Poly Coeffs:  [12, 2, 0, 12, 0, 11, 0, 6]
```

Finally, let's **decode** the new polynomial to see what happened to our messages.

```python
# E. Decode (Evaluate at original roots)
# Note: We evaluate the NEW polynomial at the OLD roots
new_msgs = [eval_poly_at_point(rotated_coeffs, r, p) for r in roots]

# F. Analysis
print("\n--- Data Movement Analysis ---")
print(f"{'Slot':<5} | {'Original (Root)':<15} | {'Moved To (Root^k)':<20} | {'New Value'}")
print("-" * 60)
for i in range(N):
    original_root = roots[i]
    target_root = shuffled_roots[i] # Where the data from this slot WENT
    
    # To understand the vector change:
    # The value at Slot i is now what WAS at the slot corresponding to target_root
    print(f"{i:<5} | {original_root:<15} | {target_root:<20} | {new_msgs[i]}")

print("\nDirect Comparison:")
print(f"Old: {msgs}")
print(f"New: {new_msgs}")
```

**Output:**
```text
--- Data Movement Analysis ---
Slot  | Original (Root) | Moved To (Root^k)    | New Value
------------------------------------------------------------
0     | 3               | 10                   | 4
1     | 5               | 6                    | 2
2     | 6               | 12                   | 6
3     | 7               | 3                    | 0
4     | 10              | 14                   | 7
5     | 11              | 5                    | 1
6     | 12              | 11                   | 5
7     | 14              | 7                    | 3

Direct Comparison:
Old: [0, 1, 2, 3, 4, 5, 6, 7]
New: [4, 2, 6, 0, 7, 1, 5, 3]
```

**What just happened?**
The vector `[0, 1, 2, 3...]` turned into `[4, 2, 6, 0...]`.
*   Look at **Slot 3** (fourth slot) (originally holding value `3` at root 7).
*   The automorphism transforms the evaluation point from $x$ to $x^3$; root 7 to root $7^3 \equiv 3 \pmod{17}$.
*   Slot 3 now looks up the value at $3 \pmod{17}$.
*   Hence, the value residing at Root 3 (which is `0`) is pulled into Slot 3.

We effectively replaced $x$ with $x^3$ to create a new polynomial (whose coefficients are `rotated_coeffs`). As a result, looking up the value at root $r_i$ in this new polynomial is mathematically equivalent to pulling the value from $r_i^3$ in the original encoding polynomial.

As can be seen in the example above, the data did shuffle, but the result looks chaotic. This is because we ordered our roots ascendingly (`3, 5, 6...`). To make this a clean "Shift Left" or "Shift Right", we need to be smarter about the order. We need the **Cyclic Group**.

> **Important Note: Fixed Slots vs. Moving Data**
> Before we rotate, remember: The **Slots** (indices 0, 1, 2...) and their assigned **Roots** never change positions. Slot 0 is always Root 3, as we specified in the setup.
> When we say "Rotation", we are not rearranging the slots themselves. We are changing the polynomial so that the **data** flows from one slot to another.

***

## Controlling the Chaos: From Galois to Cyclic Groups

In the previous example, we saw that applying a Galois automorphism ($x \to x^3$) moved our data around, but the result felt chaotic.

In most real-world FHE applications (like matrix multiplication or convolution), we do not want a random shuffle. We want a structured **Rotation**:
*   **Rotate Left:** Slot $0 \leftarrow$ Slot $1 \leftarrow$ Slot $2 \dots$
*   **Rotate Right:** Slot $0 \rightarrow$ Slot $1 \rightarrow$ Slot $2 \dots$

How do we tame the "chaos" of the Galois group to perform these neat cyclic rotations? The secret lies in how we **order** our slots.

### The Mapping: The Isomorphism
First, we must formalize the relationship between the "swaps" and the math. The Galois Group is **isomorphic** to the multiplicative group of integers modulo $m$, denoted as $(\mathbb{Z}/m\mathbb{Z})^\times$.

This isomorphism tells us two things:
1.  **Valid Moves:** The only valid "shuffles" (automorphisms) correspond to integers $k$ where $\gcd(k, m) = 1$.
2.  **Group Structure:** The way these integers behave under multiplication modulo $m$ is *exactly* how the permutations behave when composed.

### The Bridge: Generators and Ordering
If we assign our roots to slots in a simple ascending order (Root 1, Root 2, Root 3...), applying a multiplication $x \to x^k$ jumps us around the list unpredictably.

However, if the group $(\mathbb{Z}/m\mathbb{Z})^\times$ is **cyclic**, there exists a special number called a **Generator** ($g$). A generator is an integer whose powers visit every valid number in the multiplicative group before repeating.

$$ \{ g^0, g^1, g^2, g^3, \dots \} \pmod m $$

**The Strategy:** Instead of ordering roots numerically, we order them **exponentially** based on the generator.
*   **Slot 0** is assigned to Root $r^{(g^0)}$
*   **Slot 1** is assigned to Root $r^{(g^1)}$
*   **Slot 2** is assigned to Root $r^{(g^2)}$
*   ...
*   **Slot $i$** is assigned to Root $r^{(g^i)}$

It is crucial to distinguish between two types of "generators" at play:
1.  **$\zeta$ (primitive $m$-th root of unity):** This is a specific field element from which all other roots of unity (the actual evaluation points where messages are stored) are derived.
2.  **Galois generators ($g \in (\mathbb{Z}/m\mathbb{Z})^\times$, e.g., 5 and -1):** These are integers that define the *structured ordering* of the roots into logical "slots". This ordering is what translates the abstract Galois automorphism $x \to x^g$ into a predictable cyclic shift of the messages across those slots.

**The Result:**
Now, consider what happens if we choose our Galois automorphism to be $k = g$ (applying $x \to x^g$).
1.  We look at **Slot $i$**, which corresponds to root $r^{(g^i)}$.
2.  The automorphism raises this to the power of $g$.
3.  The new location is $r^{(g^i) \cdot g} = r^{(g^{i+1})}$.
4.  $r^{(g^{i+1})}$ is exactly the location of **Slot $i+1$**.

By ordering the roots according to the powers of $g$, the algebraic operation "multiply by $g$" transforms into the geometric operation "shift to the next slot."

### Python Illustration: The Generator's Orbit
To see this clearly, let's look at a case where the group is perfectly cyclic. We will use **$m=17$** (a prime number).
The group $(\mathbb{Z}/17\mathbb{Z})^\times$ has size 16. We need a generator that visits all 16 numbers.

*(Note: In this example, we are looking at the indices (exponents) of the roots, not the roots themselves).*

```python
def get_generator_orbit(m, g):
    """
    Generates the orbit of g modulo m.
    If g is a generator, this list should contain all integers 
    coprime to m.
    """
    orbit = []
    current = 1 # g^0
    
    # We loop until we wrap back around to 1
    while True:
        orbit.append(current)
        current = (current * g) % m
        if current == 1:
            break
    return orbit

# Configuration for a Prime m
m = 17 
g = 3   # 3 is a primitive root modulo 17

print(f"--- Group Structure for m={m} ---")
print(f"Generator g = {g}")

# Calculate the cycle
orbit = get_generator_orbit(m, g)

print(f"\nOrder of Slots (Powers of g):")
print(orbit)

print(f"\nTotal Slots Covered: {len(orbit)}")
print(f"Is this a full cycle? {len(orbit) == m-1}")

# simulating the shift
print("\n--- Simulating Rotation (x -> x^g) ---")
# If we are at index `orbit[i]`, multiplying by g moves us to `orbit[i+1]`
for i in range(5): # Show first 5 steps
    current_idx = orbit[i]
    next_idx = (current_idx * g) % m
    print(f"Slot {i} (Index {current_idx:2d}) * {g} -> Index {next_idx:2d} (which is Slot {i+1})")
```

**Output:**
```text
--- Group Structure for m=17 ---
Generator g = 3

Order of Slots (Powers of g):
[1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6]

Total Slots Covered: 16
Is this a full cycle? True

--- Simulating Rotation (x -> x^g) ---
Slot 0 (Index  1) * 3 -> Index  3 (which is Slot 1)
Slot 1 (Index  3) * 3 -> Index  9 (which is Slot 2)
Slot 2 (Index  9) * 3 -> Index 27 -> 10 (which is Slot 3)
Slot 3 (Index 10) * 3 -> Index 30 -> 13 (which is Slot 4)
Slot 4 (Index 13) * 3 -> Index 39 ->  5 (which is Slot 5)
```

As you can see, multiplying by $g$ ($3$) perfectly walks us through the list of slots one by one. If we map our data to these specific indices `[1, 3, 9, 10...]` instead of `[1, 2, 3, 4...]`, a Galois automorphism becomes a perfect cyclic rotation.

**What if we want to rotate by 5 slots at once?**
We do not need to apply the "shift by 1" operation five times, which would be computationally expensive. Because of the cyclic group structure, we can calculate a specific exponent $k$ that represents a "jump of 5."

$$ k = g^{\text{shift}} \pmod m $$

In our case ($g=3$, shift=5), we get $k = 3^5 \pmod{17} = 5$. By applying the automorphism $x \to x^5$ *once*, we jump 5 slots instantly.

```python
# ... (previous code) ...

# What if we want to Rotate by 5 steps at once?
shift_amount = 5

# We compute the specific automorphism exponent k
# k = g^5 mod m
k_rot5 = pow(g, shift_amount, m) 

print(f"\n--- Simulating Rotation by {shift_amount} (x -> x^{k_rot5}) ---")
# Note: k_rot5 happens to be 5 because 3^5 = 243, and 243 % 17 = 5.
# This coincidence is specific to these numbers.

start_slot = 0
current_idx = orbit[start_slot] # Value at Slot 0 (which is 1)

# Apply the automorphism math: Index * k mod m
new_idx = (current_idx * k_rot5) % m 

print(f"Start at Slot {start_slot} (Value {current_idx})")
print(f"Apply x^{k_rot5} -> New Value is {new_idx}")

# Check where this new value lives in our ordered list
target_slot = orbit.index(new_idx)

print(f"Value {new_idx} corresponds to Slot {target_slot}")
print(f"Result: Instant jump from Slot {start_slot} to Slot {target_slot}")
```

**Output:**
```text
--- Simulating Rotation by 5 (x -> x^5) ---
Start at Slot 0 (Value 1)
Apply x^5 -> New Value is 5
Value 5 corresponds to Slot 5
Result: Instant jump from Slot 0 to Slot 5
```

By choosing the correct exponent $k$, we can "teleport" data any distance around the cycle in a single step. This is the foundation of efficient rotations in FHE.

Now, there is a catch. Most FHE libraries use $m$ as a **Power of Two** (e.g., $m=8192$) rather than a prime number to maximize computational efficiency. For powers of two, a single generator is not enough. This leads us to the "Semi-Cyclic" case which we cover next.

***

## The "Common" Case: Semi-Cyclic Rotation (Power-of-Two)

While the cyclic group examples above are elegant, they rely on $m$ being a prime number (or special integer). However, the vast majority of FHE libraries (like **OpenFHE** and **SEAL**) use $m$ as a **Power of Two** (e.g., $m=4096, 8192, 16384$).

Why? Because polynomials of degree $2^k$ allow for extremely efficient **FFT/NTT algorithms**, which are crucial for performance.

### The Conflict
This efficiency comes with a mathematical quirk. The Galois group for a power-of-two cyclotomic, $(\mathbb{Z}/2^k\mathbb{Z})^\times$, is **not cyclic**. You cannot find a single generator $g$ that visits every slot.

Take $(\mathbb{Z}/16\mathbb{Z})^\times$ for example, the elements of this group are: $$\{1,3,5,7,9,11,13,15\}$$.
Let's see if any of them is a group generator. To be a generator, the powers of that number must eventually list **all** elements in the set.

*   **Try 3:**
    *   $3^1 = 3$
    *   $3^2 = 9$
    *   $3^3 = 27 \equiv 11 \pmod 16$
    *   $3^4 = 81 \equiv 1 \pmod 16$
    *   $3^5 = 243 \equiv 3 \pmod 16$. (Cycle back).
    *   **Cycle:** $$\{1, 3, 9, 11\}$$. (Length 4). 
    *   **Result:** It missed $$\{5, 7, 13, 15\}$$. Not a generator.
  
*   **Try 5:**
    *   $5^1 = 5$
    *   $5^2 = 25 \equiv 9 \pmod{16}$
    *   $5^3 = 125 \equiv 13 \pmod{16}$
    *   $5^4 = 625 \equiv 1 \pmod{16}$
    *   $5^5 = 3125 \equiv 5 \pmod{16}$. (Cycle back). 
    *   **Cycle:** $$\{1, 5, 9, 13\}$$. (Length 4).
    *   **Result:** It missed $$\{3, 7, 11, 15\}$$. Not a generator.

*   **Try 7:**
    *   $7^1 = 7$
    *   $7^2 = 49 \equiv 1 \pmod{16}$
    *   $7^3 = 343 \equiv 7 \pmod{16}$. (Cycle back).
    *   **Cycle:** $\{1, 7\}$. (Length 2).
    *   **Result:** Clearly not a generator.

You can check the rest, but the result is the same. The maximum "cycle length" (a.k.a. *order*) of any element in this group is 4, but we have 8 elements.

**Key Takeaway:** No single element can generate the entire group. The group is **not cyclic**.

This is why we need **two** generators to navigate the group structure:
1.  We use **5** to generate the main subgroup: $$\{1, 5, 9, 13\}$$. (This forms "Row 1").
2.  We use **-1** (which is 15) to reach the missing elements. Multiplying the first group by 15 gives us:
    *   $1 \times 15 = 15$
    *   $5 \times 15 = 75 \equiv 11$
    *   $9 \times 15 = 135 \equiv 7$
    *   $13 \times 15 = 195 \equiv 3$
    *   Set: $$\{15, 11, 7, 3\}$$. (This forms "Row 2").

Together, these two generators cover the entire set $$\{1, 3, 5, 7, 9, 11, 13, 15\}$$. This illustrates that for power-of-two $m$, the structure is a product of two cycles, which we will call "Semi-Cyclic".

### The General Structure: A Product of Two Cycles

For $m = 2^k$ (where $k \geq 3$), the multiplicative group $(\mathbb{Z}/m\mathbb{Z})^\times$ is never cyclic. Instead, it decomposes into a direct product of two cyclic subgroups, which can be expressed via the isomorphism:

$$C_{m/4} \times C_2$$

This means we need **two** generators to navigate the group structure:
1.  **$g = 5$**: This generator creates the main cycles.
2.  **$h = -1$ (or $m-1$)**: This generator acts as a "swap" or reflection.

One might ask, **Why 5 but not 3?** We use 5 because it provides a cleaner split of the sub-groups, ensuring that rotations do not interfere with row placement. Explaining the precise group theory behind this choice would fill an article of its own, but in essence, 5 is chosen because its powers generate the subgroup of elements congruent to $1 \pmod 4$, which provides a clean decomposition of the group structure.

### Visualizing the Hypercube (2 x m/4 Matrix)
Because of this split, we do not think of our slots as a single long line anymore. We visualize them as a **Matrix** with 2 rows.
*   **Row 1:** Powers of 5 modulo $m$ ($\{5^0, 5^1, 5^2 \dots\}$)
*   **Row 2:** Powers of 5 modulo $m$, multiplied by $-1$, that is, ($\{-5^0, -5^1, -5^2 \dots\}$)

When we pack a vector of messages, we essentially fill Row 1 first, then Row 2.

### Semi-Cyclic Rotation
This structure changes how rotations work:
*   **Applying $x \to x^5$:** This rotates the data **within the rows**.
    *   Row 1 rotates left.
    *   Row 2 rotates left.
    *   *Crucially:* Data does **not** cross from Row 1 to Row 2 automatically. This is why it is called "Semi-Cyclic."
*   **Applying $x \to x^{-1}$:** This swaps data **between the rows**.

### Python Illustration: The Two-Row Matrix
Let's visualize this with $m=32$. This gives us $N=16$ slots.
We will arrange the valid indices (the Galois group) into our 2-row structure and see how the generators affect them.

```python
def generate_row(m, start_val, g, length):
    row = []
    current = start_val
    for _ in range(length):
        row.append(current)
        current = (current * g) % m
    return row

# Parameters
m = 32
N = m // 2  # 16 slots total
row_len = N // 2 # 8 columns

# 1. Build the Matrix Structure
# Generator g=5 creates the cycle for the rows
g = 5 
# Generator h=-1 (31) jumps between rows

# Row 1 starts at 1
row1 = generate_row(m, 1, g, row_len)
# Row 2 starts at -1 (which is m-1)
row2 = generate_row(m, m-1, g, row_len)

print(f"--- Power-of-Two Structure (m={m}) ---")
print(f"Row 1 (Powers of 5):      {row1}")
print(f"Row 2 (Powers of 5 * -1): {row2}")

# 2. Simulate Rotation by 1 (x -> x^5)
print("\n--- Apply x -> x^5 (Shift within Rows) ---")
# If we multiply every index by 5, where does it go?
rotated_row1 = [(x * 5) % m for x in row1]
rotated_row2 = [(x * 5) % m for x in row2]

print(f"Row 1 moved to:           {rotated_row1}")
print(f"Row 2 moved to:           {rotated_row2}")

# Notice: The value 1 (Row1, Col0) moved to 5 (Row1, Col1).
# The rows rotated independently. Data did not cross rows.

# 3. Simulate Swap (x -> x^-1)
print("\n--- Apply x -> x^-1 (Swap Rows) ---")
# Multiply by -1 (m-1)
swapped_row1 = [(x * (m-1)) % m for x in row1]

print(f"Row 1 moved to:           {swapped_row1}")
print(f"Does this match Row 2?    {swapped_row1 == row2}")
```

**Output:**
```text
--- Power-of-Two Structure (m=32) ---
Row 1 (Powers of 5):      [1, 5, 25, 29, 17, 21, 9, 13]
Row 2 (Powers of 5 * -1): [31, 27, 7, 3, 15, 11, 23, 19]

--- Apply x -> x^5 (Shift within Rows) ---
Row 1 moved to:           [5, 25, 29, 17, 21, 9, 13, 1]
Row 2 moved to:           [27, 7, 3, 15, 11, 23, 19, 31]

--- Apply x -> x^-1 (Swap Rows) ---
Row 1 moved to:           [31, 27, 7, 3, 15, 11, 23, 19]
Does this match Row 2?    True
```

Notice that applying $x^5$ caused the lists to cycle perfectly, but Row 1 stayed in Row 1 and Row 2 stayed in Row 2, whereas applying the action $x^{-1}$ swapped the rows.

***

## Comprehensive Python Walkthrough

Let's put everything together in a single simulation. We will mock up a standard BGV/BFV rotation scenario using a **Power-of-Two** modulus ($m=16$).

Crucially, instead of sorting our roots numerically, we will order them according to the **Two-Generator Structure** (Rows) used by actual FHE libraries. This will reveal the clean "Semi-Cyclic" behavior where data rotates within the two rows of the packed matrix.

*(Note: The helper functions `inverse_crt_interpolate`, `apply_automorphism`, etc., are defined in the **[associated GitHub repository](https://github.com/caesaretos/bgv_bfv_slots_packing)** for this article.)*

```python
# 1. SETUP PARAMETERS
m = 16       # Cyclotomic index (Power of Two)
p = 17       # Modulus
N = m // 2   # Number of slots
zeta = 3     # A primitive 16-th root of unity modulo 17

print(f"--- SETUP: m={m} (Power of 2), p={p} ---")

# 2. ORGANIZE ROOTS (The "Batch" Ordering)
# Standard FHE ordering uses generators g=5 and h=-1
# Row 1: Powers of 5 [1, 5, 9, 13]
row1_exponents = [pow(5, i, m) for i in range(N // 2)]
# Row 2: Powers of 5 * -1 [15, 11, 7, 3]
row2_exponents = [(e * -1) % m for e in row1_exponents]

ordered_exponents = row1_exponents + row2_exponents
roots = [pow(zeta, e, p) for e in ordered_exponents]

print(f"Row 1 Indices: {row1_exponents}")
print(f"Row 2 Indices: {row2_exponents}")
print(f"Ordered Roots: {roots}")

# 3. ENCODE DATA
# We pack [0..3] into Row 1, and [4..7] into Row 2
msgs = [0, 1, 2, 3, 4, 5, 6, 7]
print(f"\nOriginal Data:      {msgs}")

# Pack into polynomial
points_msgs = list(zip(roots, msgs))
poly_sym = inverse_crt_interpolate(points_msgs, symbols('x'), p)
coeffs = [int(poly_sym.coeff_monomial(symbols('x')**i)) % p for i in range(N)]

# 4. PERFORM ROTATE LEFT (x -> x^5)
# This should rotate 'Left' within each row
k_rot = 5
left_rotated_coeffs = apply_automorphism(coeffs, k_rot, N, p)

# Decode
decoded_rot = [eval_poly_at_point(left_rotated_coeffs, r, p) for r in roots]
print(f"\n--- After Rotate Left (x -> x^5) ---")
print(f"Decoded Data:       {decoded_rot}")
print("Observation:        [0,1,2,3] -> [1,2,3,0] (Row 1 Rotated Left by 1)")
print("                    [4,5,6,7] -> [5,6,7,4] (Row 2 Rotated Left by 1)")

# 5. PERFORM ROW SWAP (x -> x^-1)
# This should swap Row 1 and Row 2 completely
k_swap = pow(-1, 1, m) # m-1 = 15
swapped_coeffs = apply_automorphism(coeffs, k_swap, N, p)

# Decode
decoded_swap = [eval_poly_at_point(swapped_coeffs, r, p) for r in roots]
print(f"\n--- After Swap (x -> x^-1) ---")
print(f"Decoded Data:       {decoded_swap}")
print("Observation:        Rows swapped positions completely.")

# 6. PERFORM ROTATE RIGHT (x -> x^(5^-1))
# To rotate right, we use the inverse of the generator modulo m
# g = 5. We need g_inv such that 5 * g_inv = 1 mod 16.
# 5 * 13 = 65 = (4*16) + 1. So g_inv = 13.
k_right = pow(5, -1, m) # 13
right_rot_coeffs = apply_automorphism(coeffs, k_right, N, p)

# Decode
decoded_right = [eval_poly_at_point(right_rot_coeffs, r, p) for r in roots]
print(f"\n--- After Rotate Right (x -> x^{k_right}) ---")
print(f"Decoded Data:       {decoded_right}")
print("Observation:        [0,1,2,3] -> [3,0,1,2] (Row 1 Rotated Right by 1)")
print("                    [4,5,6,7] -> [7,4,5,6] (Row 2 Rotated Right by 1)")
```

**Output:**
```text
--- SETUP: m=16 (Power of 2), p=17 ---
Row 1 Indices: [1, 5, 9, 13]
Row 2 Indices: [15, 11, 7, 3]
Ordered Roots: [3, 5, 14, 12, 6, 7, 11, 10]

Original Data:      [0, 1, 2, 3, 4, 5, 6, 7]

--- After Rotate Left (x -> x^5) ---
Decoded Data:       [1, 2, 3, 0, 5, 6, 7, 4]
Observation:        [0,1,2,3] -> [1,2,3,0] (Row 1 Rotated Left by 1)
                    [4,5,6,7] -> [5,6,7,4] (Row 2 Rotated Left by 1)

--- After Swap (x -> x^-1) ---
Decoded Data:       [4, 5, 6, 7, 0, 1, 2, 3]
Observation:        Rows swapped positions completely.

--- After Rotate Right (x -> x^13) ---
Decoded Data:       [3, 0, 1, 2, 7, 4, 5, 6]
Observation:        [0,1,2,3] -> [3,0,1,2] (Row 1 Rotated Right by 1)
                    [4,5,6,7] -> [7,4,5,6] (Row 2 Rotated Right by 1)
```

This example illustrates that by mapping the beauty of Galois automorphisms onto a carefully structured "Hypercube" of roots, we can achieve precise, programmable data movement in the encoding polynomial in BFV and BGV FHE schemes.

***

## Technical Reality: The Cost of Rotation in FHE

Before wrapping up, we must address the engineering reality. We have described rotation as a mathematical variable substitution, but does that make it slow in FHE? The answer is nuanced: the polynomial math is incredibly fast, but the cryptographic cleanup is very expensive.

### The Easy Part: Polynomial Permutation?
First, one might wonder: *Is applying automorphism, raising a polynomial to a potentially huge power like $x^{12345}$, computationally expensive?*

Surprisingly, the answer is **no**. In the context of ring arithmetic, this operation is incredibly cheap.

In FHE libraries, a polynomial is stored simply as an array of coefficients (or values): `[a0, a1, a2, ...]`. When we apply the map $x \to x^k$, we are not performing complex multiplications. We are simply **permuting** this array.
*   The term $a_1 x^1$ becomes $a_1 x^k$.
*   The term $a_2 x^2$ becomes $a_2 x^{2k}$.

The coefficient $a_i$ is simply moved to the new index $(i \cdot k) \pmod N$. If the index wraps around the modulus, we might flip the sign (because $x^N = -1$), but the core operation is just **moving memory addresses**. There is no "heavy math" involved in the polynomial transformation itself, it is merely a shuffle. 

This holds true regardless of whether the polynomial is represented in the **coefficient domain** or the **evaluation (NTT) domain**. In either case, the Galois automorphism corresponds to a specific, linear permutation of the array elements.

### Why Rotations are Expensive in FHE?
If the polynomial shuffle is so fast, why is Rotation often cited as one of the slowest operations in FHE?

The answer lies not in the *data*, but in the *secret key*. While we can easily shuffle the ciphertext polynomials, doing so fundamentally alters the mathematical relationship between the ciphertext and the secret key. We have shuffled the data as we like within the ciphertext, but the modified ciphertext becomes incompatible with the original secret key, rendering decryption impossible without a heavy maintenance operation, known as key switching.

### The Hard Part: Key Mismatch
To understand this, recall that a BGV/BFV ciphertext is a pair of polynomials $(c_0, c_1)$ that hides a message $m(x)$ relative to a secret key $s(x)$:

$$ c_0(x) + c_1(x) \cdot s(x) \approx m(x) $$

When we apply the automorphism $x \to x^k$ to the ciphertext components $c_0$ and $c_1$, the equation transforms into:

$$ c_0(x^k) + c_1(x^k) \cdot s(x^k) \approx m(x^k) $$

Do you see the problem? The decryption logic now requires the **permuted secret key**, $s(x^k)$.
However, the user only holds the original secret key $s(x)$. They cannot decrypt this rotated ciphertext!

To fix this, the server must perform a heavy maintenance operation called **Key Switching**. This requires a special public "Rotation Key" (a.k.a. Galois Key) given by the client. The server uses this key to mathematically transform the encryption from being under $s(x^k)$ back to being under $s(x)$.

**Key Takeaway:** The rotation of the data is instant; the heavy cost comes from transforming the secret key back to a usable state.

### The Implication
1.  **Storage:** A unique Rotation Key is required for every specific shift amount you want to execute in a single step (e.g., separate keys for Rotate-1, Rotate-2, Rotate-4). Storing keys for every possible shift is memory-prohibitive (often gigabytes). While one can save memory by composing smaller rotations to reach a target rotation amount (e.g., applying "Rotate 1" five times), this drastically degrades performance, as each step incurs a costly key-switching operation.
2.  **Noise:** The Key Switching process involves complex multiplications that add "noise" to the ciphertext. If you rotate too many times, the noise budget runs out, and the data becomes corrupted (decryption fails).
3.  **Speed:** Because of the key switching step, Rotation is often **orders of magnitude slower** than a simple homomorphic addition.

Therefore, while Galois rotations allow us to move data freely, they should be used sparingly in FHE applications!

***

## Key Takeaways

In this article, we bridged the gap between static data packing and dynamic data movement. Here is what you should remember:

*  **Galois Theory is the engine:** The ability to move data comes from the "indistinguishability" of roots. The automorphism map is the mathematical lever that physically moves data from one slot to another.
*  **Generators provide control:** An arbitrary automorphism index creates a random shuffle. By ordering our slots based on powers of a generator, we turn the automorphism into a clean **Cyclic Shift**.
*  **Power-of-Two means Semi-Cyclic:** In efficient instantiations of FHE schemes (where $m=2^k$), we cannot rotate the whole vector in one loop. The data is split into two batches (rows). A standard rotation shifts data *within* these rows.
*  **Automorphism is Cheap but Rotation is Expensive:** Automorphism is just data shuffling in memory, but every rotation changes the underlying secret key of the ciphertext ($s(x) \to s(x^k)$). To fix this, we must perform **Key Switching**, which adds noise and computational cost.

***

## References & Further Reading

If you want to dive deeper into the mathematics and proofs behind these concepts, here are the foundational papers that established the standard for SIMD and Ring-based BGV and BFV schemes:

1.  **The SIMD Breakthrough (Packing):** <a id="ref-simd"></a>
    N. P. Smart and F. Vercauteren, *"Fully Homomorphic SIMD Operations"*. Designs, Codes and Cryptography, 2014.
    *   [Read on IACR ePrint (2011/133)](https://eprint.iacr.org/2011/133)

2. **More on SIMD Packing and Automorphisms:** <a id="ref-simd-auto"></a>
   C. Gentry, S. Halevi and N. P. Smart, *"Fully homomorphic encryption with polylog overhead"*. EUROCRYPT 2012.
   *   [Read on IACR ePrint (2011/566)](https://eprint.iacr.org/2011/566.pdf)

3.  **The BGV Scheme:** <a id="ref-bgv"></a>
    Z. Brakerski, C. Gentry, and V. Vaikuntanathan, *"Fully Homomorphic Encryption without Bootstrapping"*. TOCT 2014.
    *   [Read on IACR ePrint (2011/277)](https://eprint.iacr.org/2011/277)

4.  **The BFV Scheme:** <a id="ref-bfv"></a>
    J. Fan and F. Vercauteren, *"Somewhat Practical Fully Homomorphic Encryption"*. IACR ePrint 2012.
    *   [Read on IACR ePrint (2012/144)](https://eprint.iacr.org/2012/144)

***

## Suggested Citation

If you found this article useful and wish to cite it in your work, we suggest:

```text
Ahmad Al Badawi, The Math of Slot Rotation in BGV and BFV, 2025, https://ahmadalbadawi.com/posts/2025/12/slot-rotation-bgv-bfv-fhe/
```

Or in BibTeX:

```bibtex
@misc{bgvbfvslotrot2025,
  author = {Al Badawi, Ahmad},
  title = {The Math of Slot Rotation in BGV and BFV},
  year = {2025},
  month = {December},
  howpublished = {\url{https://ahmadalbadawi.com/posts/2025/12/slot-rotation-bgv-bfv-fhe/}},
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
