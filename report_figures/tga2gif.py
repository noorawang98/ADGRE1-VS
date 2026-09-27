#!/usr/bin/env python3
"""Convert rendered VMD TGA frames into GIF/mp4-ready sequences.

  ~/anaconda3/envs/scvi-pypi/bin/python figures/anim/tga2gif.py            # all frame dirs
  ... tga2gif.py frames_spin figures/anim/spin.gif 15                      # one dir
"""
import glob, os, sys
import numpy as np
import imageio.v2 as imageio
from PIL import Image

def convert(frame_dir, out_gif, fps=12, max_w=700):
    fs = sorted(glob.glob(os.path.join(frame_dir, "*.tga")))
    if not fs:
        print("  no frames in", frame_dir); return None
    frames = []
    for f in fs:
        im = Image.open(f).convert("RGB")
        if im.width > max_w:
            im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
        frames.append(np.asarray(im))
    imageio.mimsave(out_gif, frames, duration=1.0 / fps, loop=0)
    print("  %s -> %s (%d frames, %.1f KB)" % (frame_dir, out_gif, len(frames),
                                                os.path.getsize(out_gif) / 1024))
    return out_gif

if __name__ == "__main__":
    if len(sys.argv) >= 3:
        convert(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 12)
    else:
        base = os.path.dirname(os.path.abspath(__file__))
        for d, name in (("frames_spin", "anim_binding_site_spin.gif"),
                        ("frames_morph", "anim_five_candidates.gif")):
            p = os.path.join(base, d)
            if os.path.isdir(p):
                convert(p, os.path.join(base, name))
