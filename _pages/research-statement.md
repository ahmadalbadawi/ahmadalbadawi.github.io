---
layout: archive
title: "Research Statement"
permalink: /research/statement/
author_profile: true
description: "Research statement by Dr. Ahmad Al Badawi: making fully homomorphic encryption fast, ergonomic, and infrastructure-ready for privacy-preserving AI at production scale."
keywords: "research statement, fully homomorphic encryption, FHE, privacy-preserving AI, lattice cryptography, OpenFHE, DARPA DPRIVE, encrypted LLM inference"
---

<p class="intro">
  The agenda behind fifteen years of work on fully homomorphic encryption, privacy-preserving machine learning, and the silicon that makes them practical.
</p>

## Private computation by default.

Most compute today runs in the clear. Training data sits in a cluster, inference servers see raw inputs, and analytics pipelines need column-level access. As models grow and inference moves into regulated industries such as healthcare, finance, and defense, that assumption gets harder to satisfy.

Fully homomorphic encryption is the most general tool we have for keeping data encrypted for its whole lifecycle, computation included. Many use cases have already shown that FHE works. What is still open is whether it can be made fast enough.

I have worked on that gap from three directions: algorithms, systems, and applications.

## 1. Algorithms

I implemented the [RNS variants of the BFV scheme](https://ieeexplore.ieee.org/abstract/document/8657794) (IEEE TETC). They ran roughly two orders of magnitude faster than CPU baselines, and they are now standard building blocks in every modern FHE library. The [first homomorphic CNN on GPUs (HCNN)](https://ieeexplore.ieee.org/abstract/document/9160866) (IEEE TETC) made encrypted deep learning something you could benchmark and measure.

[CareNets](https://priml-workshop.github.io/priml2019/papers/PriML2019_paper_28.pdf) (PriML workshop at NeurIPS 2019) took that result to high-resolution images. Its compact packing scheme fits CNN inputs, weights, and activations into HE ciphertexts. On encrypted 96x96 and 256x256 retinal images it gave over 32x speedup, 45x better memory efficiency, and a 5851x reduction in transferred message size, while staying within 3% of plaintext accuracy. [PrivFT](https://ieeexplore.ieee.org/abstract/document/9296754) (IEEE Access) carried the work over to text, training and serving private classification on encrypted documents.

## 2. Systems

Faster algorithms are only part of what deployment needs. My systems work began with the [first single-GPU CUDA implementation of FV/BFV](https://tches.iacr.org/index.php/TCHES/article/view/875) (IACR TCHES 2018), the foundational GPU-FHE result that set the baseline for every accelerated FHE implementation since. The [multi-GPU extension](https://ieeexplore.ieee.org/abstract/document/9185077) (IEEE TPDS) scaled FHE workloads across GPU clusters, so encrypted training and inference could run as parallel workloads instead of being bottlenecked on a CPU.

I co-authored [PALISADE](https://gitlab.com/palisade) and its successor [OpenFHE](https://dl.acm.org/doi/abs/10.1145/3560827.3563379), the open-source libraries that most FHE research and product work depends on today. As technical lead on DARPA DPRIVE, I took a $15M effort from concept to a 12 nm, 1 GHz ASIC design for homomorphic ML, with a custom ISA built around FHE. [TREBUCHET](https://arxiv.org/abs/2304.05237) (GOMACTech 2025) documents the architecture and methodology. Across these projects the systems work spans compilers, GPU clusters, and silicon.

## 3. Applications

The applications work is aimed at what FHE can actually deploy.

In healthcare, cross-institution oncology analysis under multiparty HE ([PNAS](https://www.pnas.org/doi/10.1073/pnas.2304415120)) showed that federated learning on real cancer-center data can be done end-to-end. CKKS-based [private pathological assessment](https://link.springer.com/article/10.1186/s13040-024-00379-9) (Springer BioData Mining, 2024) applied the same machinery to encrypted pathology classification. SVM inference on encrypted patient data, paired with a compact feature-extraction pipeline, runs in seconds at 128-bit security and matches the accuracy of plaintext baselines.

For financial-sector analytics, [FHSVM](https://link.springer.com/article/10.1007/s00521-022-07202-8) (*Neural Computing and Applications*, 2022) ran homomorphic SVM inference for anti-money-laundering classification on encrypted Bitcoin-transaction datasets. Using novel CKKS packing strategies and a parallel implementation, it reached roughly 1.25 s prediction latency on multi-core CPUs at 128-bit security, with zero accuracy loss compared to the plaintext model.

For LLMs, I introduced [POLARIS](https://eprint.iacr.org/2026/935), an open-source, model-preserving reference framework for CKKS-based private LLM inference. It runs encrypted BERT-Tiny, BERT-Mini, BERT-Small, BERT-Medium, and BERT-Base with GPU acceleration and is evaluated on standard, unmodified architectures, with no retraining and no activation-function substitution. POLARIS is meant as a proof-of-concept and a shared baseline for the community. It comes with the Private LLM Card System (PLCS), a standard way to report framework configurations and results across the field.

Cancer-center collaborations, federated learning on encrypted gradients, intrusion detection on encrypted telemetry, and private LLM inference on encrypted prompt embeddings together show what production-grade private computation looks like.

## The next decade

The frontier is private inference for large generative models. CKKS-based frameworks have now run end-to-end encrypted inference on LLMs of up to 8 billion parameters. My recent systematization of knowledge ([SoK on private LLM inference under approximate HE](https://eprint.iacr.org/2026/935)) finds a runtime gap of roughly four orders of magnitude between encrypted and plaintext inference, and identifies it as the main barrier to practical use. CKKS-based inference is algorithmically feasible today on standard, unmodified LLMs from BERT-Tiny up to Llama-3-8B, but it will not be operationally practical for human-facing applications until that gap narrows.

Closing the gap needs progress on all three pillars at the same time:

- Algorithms: FHE-friendly approximations of transformer building blocks (softmax, GELU, layer norm), packing layouts for linear and non-linear blocks, and bootstrapping schedules tuned to attention and MLP instead of generic worst-case loops.
- Systems: compilers, hybrid execution stacks, GPU acceleration, and silicon that turn algorithmic gains into wall-clock latency a user will tolerate.
- Applications: regulated, latency-tolerant use cases (healthcare diagnostics, federated analytics, multi-party machine learning, secure inference for defense) where the privacy guarantee is worth a real latency budget today, and where deployment exercises the rest of the stack.

A 10x improvement on each of the three pillars compounds to a 1000x reduction. That closes most of the four-orders-of-magnitude gap and leaves encrypted inference roughly 10x slower than plaintext. Latency-tolerant production workloads can absorb that overhead easily: overnight analytics, federated training rounds, regulated medical and financial diagnostics, batch scoring, and asynchronous LLM serving. At that point the privacy guarantee starts to outweigh the cost. Each of the three 10x gains is within reach on its own. The open question is whether the field can deliver them together, instead of each pillar shipping its own 10x in isolation.

I work across all three pillars, and I am looking for collaborators (labs, startups, program committees) who want to make encrypted AI practical.

<p class="section-more">
  <a class="btn btn--primary" href="/contact/">Get in touch</a>
  <a class="btn btn--ghost" href="/research/">See current research</a>
  <a class="btn btn--ghost" href="/research/impact/">See research impact</a>
</p>
