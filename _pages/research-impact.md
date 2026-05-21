---
layout: archive
title: "Research Impact"
permalink: /research/impact/
author_profile: true
description: "Real-world adoption of Dr. Ahmad Al Badawi's research: OpenFHE leading the FHE library market, RNS-BFV as standard library infrastructure, the DPRIVE FHE ASIC, and papers that have crossed 100 citations."
---

<p class="intro">
  Citations measure attention; this page tracks what the work has actually been <em>used</em> for. Libraries shipped and adopted across the field, hardware taped out, patient data analyzed across hospitals, and the papers that have crossed milestone citation counts.
</p>

## OpenFHE leads the FHE library market

<article class="feature-card">
  <header class="feature-card__meta">
    <span class="badge badge--accent">External signal</span>
    <span class="badge">April 2026</span>
  </header>
  <h3 class="feature-card__title">51% adoption. #1 fully homomorphic encryption library.</h3>
  <div class="feature-card__body">
    <p>Lattica's <a href="https://www.lattica.ai/news/2nd-fhe-landscape-survey" target="_blank" rel="noopener noreferrer">2nd FHE Landscape Survey</a> (16 April 2026) of cryptographers, engineers, and researchers ranks <strong>OpenFHE</strong> as the most widely used fully homomorphic encryption library, ahead of TFHE-rs (33%), Lattigo (21%), Concrete (18%), HeaAN (13%), HElayers (10%), SEAL/TenSEAL (10%), HELib (8%), and Poulpy (5%).</p>
    <p><a href="https://github.com/openfheorg/openfhe-development" target="_blank" rel="noopener noreferrer">OpenFHE</a> is the open-source FHE library I co-authored, the successor to <a href="https://gitlab.com/palisade" target="_blank" rel="noopener noreferrer">PALISADE</a>. The signal is strengthened by the fact that Lattica itself builds a competing commercial FHE library (listed in the same survey under "proprietary library, 15%"), so this is an independent measurement rather than a self-report.</p>
  </div>
</article>

## Footprint

<ul class="card-grid">
  <li>
    <a class="card" href="https://github.com/openfheorg/openfhe-development" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">PALISADE to OpenFHE - a decade of open-source FHE</h3>
      <div class="card__meta">
        <span>Co-author</span>
        <span>Used across academia and industry</span>
      </div>
      <p class="card__excerpt">Co-authored both PALISADE and its successor OpenFHE. Together they represent over a decade of open-source FHE infrastructure development, now standard tooling across the FHE research community and a growing number of commercial deployments.</p>
    </a>
  </li>

  <li>
    <a class="card" href="https://doi.org/10.1109/TETC.2019.2902799" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">RNS-BFV (BEHZ + HPS) - standard primitives</h3>
      <div class="card__meta">
        <span>IEEE TETC 2019</span>
        <span>Adopted by multiple FHE libraries</span>
      </div>
      <p class="card__excerpt">The BEHZ and HPS RNS variants of BFV introduced in this paper are now baseline building blocks inside modern FHE libraries (OpenFHE, SEAL, HElib), delivering roughly two orders of magnitude speedup over the original CPU baselines.</p>
    </a>
  </li>

  <li>
    <a class="card" href="https://arxiv.org/abs/2304.05237" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">TREBUCHET / DPRIVE - first 12 nm FHE ASIC</h3>
      <div class="card__meta">
        <span>$15M DARPA</span>
        <span>1 GHz, 176 mm&sup2; floorplan</span>
      </div>
      <p class="card__excerpt">Technical lead on the DPRIVE program: a custom 12 nm ASIC for homomorphic deep learning. Delivered tape-out-ready RTL with 1 GHz timing closure, custom ISA, microcode scheduler, and node-array architecture.</p>
    </a>
  </li>
</ul>

## Real-world deployment

<article class="feature-card">
  <header class="feature-card__meta">
    <span class="badge badge--accent">PNAS</span>
    <span class="badge">2023</span>
    <span class="badge">Clinical data</span>
  </header>
  <h3 class="feature-card__title">Multiparty FHE on real oncology data</h3>
  <div class="feature-card__body">
    <p>Collaborative privacy-preserving analysis of oncological data using multiparty homomorphic encryption - cross-institution analysis of real patient records without exposing any individual hospital's data. Not a benchmark, not a synthetic dataset: an end-to-end deployment in a high-sensitivity domain, published in the Proceedings of the National Academy of Sciences.</p>
    <p><a href="https://www.pnas.org/doi/10.1073/pnas.2304415120" target="_blank" rel="noopener noreferrer">Read the paper &rarr;</a></p>
  </div>
</article>

## Papers with 100+ citations

<p class="intro">The following papers have crossed 100 citations on Google Scholar. Listed in citation-count order. Last updated, May 2026.</p>

<ul class="card-grid">
  <li>
    <a class="card" href="https://doi.org/10.1145/3560827.3563379" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">OpenFHE: Open-Source Fully Homomorphic Encryption Library</h3>
      <div class="card__meta">
        <span>WAHC</span>
        <span>2022</span>
      </div>
    </a>
  </li>

  <li>
    <a class="card" href="https://doi.org/10.1109/TETC.2020.3014636" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">Towards the AlexNet Moment for Homomorphic Encryption: HCNN, the First Homomorphic CNN on Encrypted Data with GPUs</h3>
      <div class="card__meta">
        <span>IEEE TETC</span>
        <span>2018</span>
      </div>
    </a>
  </li>

  <li>
    <a class="card" href="https://doi.org/10.1109/TETC.2019.2902799" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">Implementation and Performance Evaluation of RNS Variants of the BFV Homomorphic Encryption Scheme</h3>
      <div class="card__meta">
        <span>IEEE TETC</span>
        <span>2019</span>
      </div>
    </a>
  </li>

  <li>
    <a class="card" href="https://doi.org/10.1109/ACCESS.2020.3045465" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">PrivFT: Private and Fast Text Classification with Homomorphic Encryption</h3>
      <div class="card__meta">
        <span>IEEE Access</span>
        <span>2020</span>
      </div>
    </a>
  </li>

  <li>
    <a class="card" href="https://link.springer.com/article/10.1007/s00521-022-07015-9" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">High-Performance Intrusion Detection System for Networked UAVs via Deep Learning</h3>
      <div class="card__meta">
        <span>Neural Comput. &amp; Appl.</span>
        <span>2022</span>
      </div>
    </a>
  </li>

  <li>
    <a class="card" href="https://tches.iacr.org/index.php/TCHES/article/view/875" target="_blank" rel="noopener noreferrer">
      <h3 class="card__title">High-Performance FV Somewhat Homomorphic Encryption on GPUs: An Implementation using CUDA</h3>
      <div class="card__meta">
        <span>IACR TCHES</span>
        <span>2018</span>
      </div>
    </a>
  </li>
</ul>

<p class="section-more"><a href="/publications/">All publications &rarr;</a></p>
