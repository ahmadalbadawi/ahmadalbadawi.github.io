# ahmadalbadawi.com

Personal website of Dr. Ahmad Al Badawi - Senior Scientist, Cryptography.

Custom Jekyll theme built in place, deployed via GitHub Pages on the custom
domain [ahmadalbadawi.com](https://ahmadalbadawi.com). Originally forked from
Academic Pages / Minimal Mistakes, then rewritten with modern CSS, vanilla JS,
WebP imagery, dark mode, and accessibility features.

## Develop locally

```bash
bundle install
bundle exec jekyll serve --config _config.yml,_config.dev.yml
```

Open http://localhost:4000 in a browser. Auto-reloads on file changes.

## Add content

| Type        | Where                                       |
|-------------|---------------------------------------------|
| Blog post   | `_posts/YYYY-MM-DD-slug.md`                 |
| Talk        | `_talks/YYYY-MM-DD-slug.md`                 |
| Teaching    | `_teaching/YYYY-AB-cd.md`                   |
| Portfolio   | `_portfolio/slug.md`                        |
| Static page | `_pages/slug.md` (add to `_data/navigation.yml` for menu) |

Each file needs YAML front matter with at least `title:` and `date:`. For pages
with math, add `mathjax: true` to load MathJax 3 only on that page.

## Optimize new images

The site uses `<picture>` with WebP + JPG fallback, plus a 480px square crop
for avatars. Whenever you add an image to `/images/`, regenerate the variants:

```bash
script/optimize-images.sh
```

Requires ImageMagick (`apt install imagemagick` on Debian/Ubuntu).

## Deploy

Push to `master`. GitHub Pages builds and deploys via the `github-pages` gem
(see `Gemfile`). No GitHub Actions workflow required.

The `CNAME` file sets the custom domain to `ahmadalbadawi.com`. DNS A/AAAA
records point at GitHub Pages.

## Theme & accessibility notes

- **Dark mode:** auto-follows `prefers-color-scheme`; toggle in the header
  overrides and persists in `localStorage` (cross-tab synced).
- **No jQuery / FontAwesome.** Two ~1 KB vanilla JS files, inline SVG icons,
  system font stacks. CSS is ~20 KB gzipped.
- **Skip link, focus-visible styles, `aria-*` labels** on all icon-only
  controls, `prefers-reduced-motion` respected.
- **SEO:** custom `_includes/seo.html` emits OG/Twitter cards (with image
  dimensions), Person + Article JSON-LD, canonical URLs, sitemap, atom feed.

## Configuration knobs

In `_config.yml`:

- `site.author.*` drives the sidebar profile card and JSON-LD `sameAs` links.
- `site.og_image*` controls the default social-share card (`/images/og-card.jpg`).
- `analytics.provider: "google-gtag"` + `analytics.google.tracking_id: G-XXX`
  enables Google Analytics 4. Leave `false` to skip analytics.

## Project layout

```
_config.yml              # site config
_config.dev.yml          # localhost overrides
_data/navigation.yml     # top nav menu items
_includes/               # liquid partials
_layouts/                # page layouts (default, single, archive, talk)
_pages/                  # static pages (about, education, …)
_posts/                  # dated blog entries
_portfolio/, _talks/, _teaching/   # collection content
_sass/                   # source SCSS partials (compiled to assets/css/main.css)
assets/css/main.scss     # one-liner entry that imports _sass/main
assets/js/               # theme.js, nav.js, copy-code.js
images/                  # site imagery (+ WebP variants)
script/optimize-images.sh
```
