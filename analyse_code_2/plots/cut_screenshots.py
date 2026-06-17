#!/usr/bin/env python3
import sys
from pathlib import Path

from PIL import Image  # pip install pillow

CROP_SIZE = 400  # pixels


def center_crop_500(img: Image.Image) -> Image.Image:
    """Return a 500x500 center crop. Raises ValueError if image is too small."""
    w, h = img.size

    if w < CROP_SIZE or h < CROP_SIZE:
        raise ValueError(f"Image is too small for a {CROP_SIZE}x{CROP_SIZE} crop (size: {w}x{h})")

    left = (w - CROP_SIZE) // 2
    top = (h - CROP_SIZE) // 2
    right = left + CROP_SIZE
    bottom = top + CROP_SIZE

    return img.crop((left, top, right, bottom))


def process_image(path: Path) -> None:
    try:
        with Image.open(path) as img:
            img = img.convert("RGB")  # ensure JPEG-compatible
            cropped = center_crop_500(img)

            out_path = path.with_name(f"{path.stem}_center_cut.jpg")
            cropped.save(out_path, format="JPEG", quality=95)
            print(f"Saved center crop: {out_path}")
    except Exception as e:
        print(f"Error processing {path}: {e}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python center_cut.py image1.png image2.jpg ...")
        sys.exit(1)

    for arg in sys.argv[1:]:
        path = Path(arg)
        if not path.is_file():
            print(f"Not a file or does not exist: {path}")
            continue
        process_image(path)


if __name__ == "__main__":
    main()