#!/usr/bin/env python3
import sys
from pathlib import Path

from PIL import Image  # pip install pillow

CROP_SIZE_X = 400  # pixels
CROP_SIZE_Y = 400  # pixels


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


def crop_with_margins(
    img: Image.Image,
    cut_left: int,
    cut_top: int,
    cut_right: int,
    cut_bottom: int
) -> Image.Image:
    """
    Crop the image by cutting a given number of pixels from each side.

    :param img: Input PIL image.
    :param cut_left: Pixels to remove from the left side.
    :param cut_top: Pixels to remove from the top side.
    :param cut_right: Pixels to remove from the right side.
    :param cut_bottom: Pixels to remove from the bottom side.
    :return: Cropped PIL image.
    :raises ValueError: If the crop would result in non-positive dimensions.
    """
    w, h = img.size

    new_width = w - cut_left - cut_right
    new_height = h - cut_top - cut_bottom

    if new_width <= 0 or new_height <= 0:
        raise ValueError(
            f"Invalid crop: resulting size would be {new_width}x{new_height} "
            f"from original {w}x{h}."
        )

    left = cut_left
    top = cut_top
    right = w - cut_right
    bottom = h - cut_bottom

    return img.crop((left, top, right, bottom))


def process_image(path: Path, cut_left, cut_top, cut_right, cut_bottom) -> None:
    try:
        with Image.open(path) as img:
            img = img.convert("RGB")  # ensure JPEG-compatible
            cropped = crop_with_margins(img, cut_left, cut_top, cut_right, cut_bottom)

            out_path = path.with_name(f"{path.stem}_center_cut.jpg")
            cropped.save(out_path, format="JPEG", quality=95)
            print(f"Saved center crop: {out_path}")
    except Exception as e:
        print(f"Error processing {path}: {e}")


def process_path_list(path_list, cut_left, cut_top, cut_right, cut_bottom):
    for item in path_list:
        path = Path(item)
        if not path.is_file():
            print(f"Not a file or does not exist: {path}")
            continue
        process_image(path, cut_left, cut_top, cut_right, cut_bottom)

def main():
#     if len(sys.argv) < 2:
#         print("Usage: python center_cut.py image1.png image2.jpg ...")
#         sys.exit(1)

    path_100 = ["coarsening_ovito/PVA_100_T088_tstep=0.jpg", "coarsening_ovito/PVA_100_T088_tstep=16=tc.jpg", "coarsening_ovito/PVA_100_T088_tstep=30.jpg"]
    path_1000 = ["coarsening_ovito/PVA_1000_T088_tstep=0.jpg", "coarsening_ovito/PVA_1000_T088_tstep=12=tc.jpg", "coarsening_ovito/PVA_1000_T088_tstep=30.jpg"]
    # for arg in sys.argv[1:]:
    #     path = Path(arg)

    #process_path_list(path_100, 400, 85, 650, 150)
    process_path_list(path_1000, 200, 100, 500, 150)

if __name__ == "__main__":
    main()