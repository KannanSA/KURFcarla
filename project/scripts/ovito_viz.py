# scripts/ovito_viz.py
import os
from ovito.io import import_file, export_file

# Export a snapshot image for each LAMMPS dump
def export_snapshot(dump_file, out_image, frame=0, dpi=300, size=(800,600)):
    pipeline = import_file(dump_file)
    export_file(pipeline, out_image, 'png', frame=frame, size=size, dpi=dpi)

if __name__ == '__main__':
    os.makedirs('project/results/ovito', exist_ok=True)
    for fname in os.listdir('project/out'):
        if fname.endswith('.dump'):
            in_dump = os.path.join('project/out', fname)
            out_img = os.path.join('project/results/ovito', fname.replace('.dump', '.png'))
            export_snapshot(in_dump, out_img)
