---
layout: archive
title: "Research Statement"
permalink: /research/statement/
author_profile: true
description: "Research statement by Dr. Ahmad Al Badawi: making fully homomorphic encryption fast, ergonomic, and infrastructure-ready for privacy-preserving AI at production scale."
keywords: "research statement, fully homomorphic encryption, FHE, privacy-preserving AI, lattice cryptography, OpenFHE, DARPA DPRIVE, encrypted LLM inference"
---

<p class="intro">
  The agenda behind fifteen years of work on fully homomorphic encryption, privacy-preserving machine learning, and the silicon to make them practical.
</p>

## Private computation should be a default, not a feature.

Most of today's compute runs in the clear. Training data lives in a cluster, inference servers see raw inputs, analytics pipelines need column-level access. As models scale and inference moves into regulated industries (healthcare, finance, defense), that assumption is increasingly the bottleneck.

Fully homomorphic encryption is the most general tool we have for keeping data encrypted throughout its lifecycle, including during computation. The question is no longer whether FHE works. It does. The question is whether it can be made fast, ergonomic, and infrastructure-ready.

My work has been a sustained attack on that gap, from three angles.

## 1. Algorithms

The [RNS variants of the BFV scheme](https://ieeexplore.ieee.org/abstract/document/8657794) that I implemented delivered roughly two orders of magnitude over CPU baselines and are now standard building blocks in every modern FHE library. The [first homomorphic CNN on GPUs (HCNN)](https://ieeexplore.ieee.org/abstract/document/9160866) showed that encrypted deep learning is not a thought experiment but a measurable benchmark. [CareNets](https://priml-workshop.github.io/priml2019/papers/PriML2019_paper_28.pdf) (NeurIPS 2019) pushed that result to high-resolution images: a compact packing scheme that fits CNN inputs, weights, and activations into HE ciphertexts, delivering over 32x speedup, 45x better memory efficiency, and a 5851x reduction in transferred message size on encrypted 96x96 and 256x256 retinal images, all within 3% of the plaintext accuracy. [PrivFT](https://ieeexplore.ieee.org/abstract/document/9296754) extended the line into text, training and serving private classification on encrypted documents. More recent work on [private pathological assessment](https://link.springer.com/article/10.1186/s13040-024-00379-9) pushed the same approach into medical imaging: CKKS-based SVM inference on encrypted patient data, paired with a compact feature-extraction pipeline, runs encrypted pathology classification in seconds at 128-bit security and matches the accuracy of plaintext baselines.

## 2. Systems

Algorithms alone do not close the deployment gap. My systems work started with the [first single-GPU CUDA implementation of FV/BFV](https://tches.iacr.org/index.php/TCHES/article/view/875) (IACR TCHES 2018), the foundational GPU-FHE result that set the baseline for every accelerated FHE implementation since. The [multi-GPU extension](https://ieeexplore.ieee.org/abstract/document/9185077) then scaled FHE workloads across GPU clusters, turning encrypted training and inference into a parallel workload rather than a CPU bottleneck. As a co-author of [OpenFHE](https://dl.acm.org/doi/abs/10.1145/3560827.3563379), I helped build the open-source infrastructure that most FHE research and product work today depends on. As technical lead on DARPA DPRIVE, I drove a $15M effort from concept to a 12 nm, 1 GHz ASIC design for homomorphic ML, with a custom ISA designed around FHE. The architecture and methodology are documented in [TREBUCHET](https://arxiv.org/abs/2304.05237) (GOMACTech 2025). From compilers to GPU clusters to silicon, the systems work is the rest of the answer.

## 3. Applications

From the [PNAS paper](https://www.pnas.org/doi/10.1073/pnas.2304415120) on cross-institution oncology analysis under multiparty HE to current work on private LLM inference, the goal is to push what FHE can deploy, not just what it can demonstrate. Cancer-center collaborations, federated learning on encrypted gradients, and intrusion detection on encrypted telemetry are not toy demos. They are the cases that show what production-grade private computation actually looks like.

## The next decade

The frontier is private inference for large generative models. CKKS-based frameworks have now demonstrated end-to-end encrypted inference of LLMs up to 8 billion parameters, but my recent systematization of knowledge ([SoK on private LLM inference under approximate HE](https://eprint.iacr.org/2026/935)) identifies a runtime gap of roughly **four orders of magnitude** between encrypted and plaintext inference as the primary barrier to practical use. CKKS-based inference is now algorithmically feasible on standard, unmodified LLMs from BERT-Tiny up to Llama-3-8B; it is not yet operationally practical for human-facing applications until that efficiency gap is narrowed.

Closing it is not a single-discipline problem. It requires progress on all three pillars at once:

- **Algorithms**: FHE-friendly approximations of transformer building blocks (softmax, GELU, layer norm), packing layouts for linear and non-linear blocks, and bootstrapping schedules tuned to attention and MLP rather than generic worst-case loops.
- **Systems**: compilers, hybrid execution stacks, GPU acceleration, and silicon that turn the algorithmic gains into wall-clock latency a user will tolerate.
- **Applications**: regulated, latency-tolerant use cases (healthcare diagnostics, federated analytics, multi-party machine learning, secure inference for defense) where the privacy guarantee is worth a real latency budget today, and where deployment exercises the rest of the stack.

I work across all three and I am looking for collaborators (labs, startups, program committees) who want to make encrypted-by-default AI a practical industry reality.

<p class="section-more">
  <a class="btn btn--primary" href="/contact/">Get in touch</a>
  <a class="btn btn--ghost" href="/research/">See current research</a>
</p>
