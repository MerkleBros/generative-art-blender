import sys
import os
import bpy

blend_dir = os.path.dirname(bpy.data.filepath)
if blend_dir not in sys.path:
   sys.path.append(blend_dir)

import art
import perlin
import importlib
importlib.reload(art)
importlib.reload(perlin)
art.main()
