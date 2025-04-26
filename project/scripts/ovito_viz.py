import os
from ovito.io import import_file, export_file

def export_snapshot(dump_file, out_image, frame=0, dpi=300, size=(800,600)):
    pipeline = import_file(dump_file)
    export_file(pipeline, out_image, 'png', frame=frame, size=size, dpi=dpi)

if __name__ == '__main__':
    os.makedirs('results/ovito', exist_ok=True)
    for fname in os.listdir('out/'):
        if fname.endswith('.dump'):
            in_dump = os.path.join('out/', fname)
            out_img  = os.path.join('results/ovito', fname.replace('.dump', '.png'))
            export_snapshot(in_dump, out_img)
