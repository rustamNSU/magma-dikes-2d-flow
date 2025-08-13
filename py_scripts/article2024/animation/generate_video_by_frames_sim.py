import os, sys
repository_dir = os.path.abspath(os.getcwd())
sys.path.append(repository_dir)

import glob
from moviepy import ImageSequenceClip
import imageio.v2 as imageio


simID = 161
fps = 15
output_mp4 = os.path.join(repository_dir, "images", f"simID{simID}", "animation.mp4")
output_gif = os.path.join(repository_dir, "images", f"simID{simID}", "animation.gif")
frame_dir = os.path.join(repository_dir, "images", f"simID{simID}", "frames")
frames_sorted = sorted(glob.glob(os.path.join(frame_dir, "frame_*.png")))


frames = sorted(glob.glob(os.path.join(frame_dir, "frame_*.png")))
if not frames:
    raise RuntimeError(f"Нет кадров в {frame_dir}")
print(f"Создаём видео из {len(frames)} кадров…")

clip = ImageSequenceClip(frames, fps=fps)
clip.write_videofile(output_mp4, codec='libx264', audio=False)
# clip.write_gif(output_gif, fps=fps)
print(f"Готово: {output_mp4} и {output_gif}")