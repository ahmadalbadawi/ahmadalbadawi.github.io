#!/usr/bin/env bash
# Generate optimized variants of images in /images/.
#
# For every JPG/PNG without a matching .webp, generates:
#   - <name>.webp           (same dimensions, smaller bytes)
#   - <name>-480.jpg/webp   (480x480 square cover crop, for avatar use)
#
# Run from the repo root. Requires ImageMagick (`convert`).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMG_DIR="$ROOT/images"

cd "$IMG_DIR"

shopt -s nullglob nocaseglob
for src in *.jpg *.jpeg *.png; do
  [[ "$src" == *-480.* ]] && continue
  base="${src%.*}"

  webp="$base.webp"
  if [[ ! -f "$webp" || "$src" -nt "$webp" ]]; then
    echo "  $src → $webp"
    convert "$src" -strip -quality 75 "$webp"
  fi

  avatar_jpg="$base-480.jpg"
  if [[ ! -f "$avatar_jpg" || "$src" -nt "$avatar_jpg" ]]; then
    echo "  $src → $avatar_jpg (480² crop)"
    convert "$src" -auto-orient -resize '480x480^' -gravity center -extent 480x480 \
      -strip -quality 82 -interlace JPEG "$avatar_jpg"
  fi

  avatar_webp="$base-480.webp"
  if [[ ! -f "$avatar_webp" || "$avatar_jpg" -nt "$avatar_webp" ]]; then
    echo "  $avatar_jpg → $avatar_webp"
    convert "$avatar_jpg" -strip -quality 75 "$avatar_webp"
  fi
done

echo "Done."
