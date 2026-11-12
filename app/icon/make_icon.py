#!/usr/bin/env python3
"""Generate the ΨLab app icon (macOS asset catalog) with Pillow.

Usage:  python3 app/icon/make_icon.py
Writes app/MacQC/Assets.xcassets/AppIcon.appiconset/*.png + Contents.json and
app/icon/AppIcon-1024.png (master). Design: macOS icon grid (824 pt rounded square
on a 1024 canvas), indigo-to-teal gradient, two faint orbital ellipses, white Ψ.
"""
import json
import math
import pathlib

from PIL import Image, ImageChops, ImageDraw, ImageFilter, ImageFont

ROOT = pathlib.Path(__file__).resolve().parents[2]
ICONSET = ROOT / "app/MacQC/Assets.xcassets/AppIcon.appiconset"
MASTER = ROOT / "app/icon/AppIcon-1024.png"
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf"

SS = 4                      # supersampling factor for smooth edges
N = 1024 * SS
TOP_LEFT = (27, 38, 112)    # deep indigo
BOTTOM_RIGHT = (22, 168, 170)  # teal


def superellipse_mask(size, inset, exponent=5.0):
    """Rounded-square ("squircle") mask: |x|^n + |y|^n <= 1 inside the inset box."""
    mask = Image.new("L", (size, size), 0)
    half = (size - 2 * inset) / 2
    cx = cy = size / 2
    pts = []
    for i in range(720):
        t = 2 * math.pi * i / 720
        c, s = math.cos(t), math.sin(t)
        x = cx + half * math.copysign(abs(c) ** (2 / exponent), c)
        y = cy + half * math.copysign(abs(s) ** (2 / exponent), s)
        pts.append((x, y))
    ImageDraw.Draw(mask).polygon(pts, fill=255)
    return mask


def gradient(size):
    g = Image.new("RGB", (size, size))
    px = g.load()
    for y in range(0, size, SS):
        for x in range(0, size, SS):
            t = (x + y) / (2 * size)
            col = tuple(int(a + (b - a) * t) for a, b in zip(TOP_LEFT, BOTTOM_RIGHT))
            for dy in range(SS):
                for dx in range(SS):
                    px[x + dx, y + dy] = col
    return g


def make_master():
    inset = 100 * SS                      # macOS grid: 824/1024 body
    mask = superellipse_mask(N, inset)
    body = gradient(N).convert("RGBA")

    # Orbital ellipses (faint white), drawn on their own layer and rotated.
    orbit = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    od = ImageDraw.Draw(orbit)
    w, h = 300 * SS, 120 * SS
    box = (N / 2 - w, N / 2 - h, N / 2 + w, N / 2 + h)
    od.ellipse(box, outline=(255, 255, 255, 70), width=14 * SS)
    layer = Image.alpha_composite(orbit.rotate(30, resample=Image.BICUBIC),
                                  orbit.rotate(-30, resample=Image.BICUBIC))
    body = Image.alpha_composite(body, layer)

    # Ψ glyph, centred on its ink box, with a soft shadow.
    font = ImageFont.truetype(FONT, 560 * SS)
    glyph = Image.new("L", (N, N), 0)
    gd = ImageDraw.Draw(glyph)
    l, t, r, b = gd.textbbox((0, 0), "Ψ", font=font)
    gd.text(((N - (r - l)) / 2 - l, (N - (b - t)) / 2 - t), "Ψ", font=font, fill=255)
    shadow = glyph.filter(ImageFilter.GaussianBlur(10 * SS))
    shadow_rgba = Image.new("RGBA", (N, N), (0, 0, 40, 0))
    shadow_rgba.putalpha(shadow.point(lambda v: v * 0.45))
    shadow_rgba = ImageChops.offset(shadow_rgba, 0, 10 * SS)
    body = Image.alpha_composite(body, shadow_rgba)
    white = Image.new("RGBA", (N, N), (255, 255, 255, 255))
    white.putalpha(glyph)
    body = Image.alpha_composite(body, white)

    body.putalpha(mask)
    return body.resize((1024, 1024), Image.LANCZOS)


SIZES = [(16, 1), (16, 2), (32, 1), (32, 2), (128, 1), (128, 2),
         (256, 1), (256, 2), (512, 1), (512, 2)]


def main():
    master = make_master()
    MASTER.parent.mkdir(parents=True, exist_ok=True)
    master.save(MASTER)
    ICONSET.mkdir(parents=True, exist_ok=True)
    images = []
    for pt, scale in SIZES:
        px = pt * scale
        name = f"icon_{pt}x{pt}{'@2x' if scale == 2 else ''}.png"
        master.resize((px, px), Image.LANCZOS).save(ICONSET / name)
        images.append({"idiom": "mac", "size": f"{pt}x{pt}", "scale": f"{scale}x", "filename": name})
    info = {"author": "xcode", "version": 1}
    (ICONSET / "Contents.json").write_text(json.dumps({"images": images, "info": info}, indent=2) + "\n")
    (ICONSET.parent / "Contents.json").write_text(json.dumps({"info": info}, indent=2) + "\n")
    print("wrote", MASTER, "and", len(images), "icon sizes")


if __name__ == "__main__":
    main()
