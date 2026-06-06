---
layout: archive
title: "Selected Publications & Patents"
permalink: /publications/
author_profile: true
description: "Peer-reviewed publications by Dr. Ahmad Al Badawi on fully homomorphic encryption, privacy-preserving machine learning, and lattice cryptography."
---

<p class="intro">
  A selection of journal, conference, and preprint publications.
  For the complete list, see my
  <a href="https://scholar.google.com.sg/citations?user=-EhCfyEAAAAJ&amp;hl=en" target="_blank" rel="noopener noreferrer">Google Scholar profile</a>.
</p>

{% assign papers = site.data.publications | where_exp: "p", "p.kind != 'patent'" %}
{% assign years = papers | map: "year" | uniq | sort | reverse %}

{% for year in years %}
<section class="year-section">
  <h2>{{ year }}</h2>
  <ul class="pub-list">
    {% assign year_papers = papers | where: "year", year %}
    {% for pub in year_papers %}
      {% include pub-item.html pub=pub %}
    {% endfor %}
  </ul>
</section>
{% endfor %}

<hr class="divider-wide">

## Patents

{% assign patents = site.data.publications | where: "kind", "patent" %}
{% assign patent_years = patents | map: "year" | uniq | sort | reverse %}
<ul class="pub-list">
  {% for y in patent_years %}
    {% assign year_patents = patents | where: "year", y %}
    {% for pub in year_patents %}
      {% include pub-item.html pub=pub %}
    {% endfor %}
  {% endfor %}
</ul>
