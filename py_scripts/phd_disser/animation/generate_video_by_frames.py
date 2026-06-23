import os, sys
repository_dir = os.path.abspath(os.getcwd())
sys.path.append(repository_dir)

import glob
from moviepy import ImageSequenceClip
import imageio.v2 as imageio


fps = 15
# dir_name = "temperature_frames_117_110"
dir_name = "beta_frames_117_110"
frame_dir = os.path.join(repository_dir, "images", dir_name)
output_mp4 = os.path.join(repository_dir, "images", dir_name + ".mp4")
frames = sorted(glob.glob(os.path.join(frame_dir, "frame_*.png")))
if not frames:
    raise RuntimeError(f"Нет кадров в {frame_dir}")
print(f"Создаём видео из {len(frames)} кадров…")

clip = ImageSequenceClip(frames, fps=fps)
clip.write_videofile(output_mp4, codec='libx264', audio=False)
# clip.write_gif(output_gif, fps=fps)
print(f"Готово: {output_mp4}")