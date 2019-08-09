import yaml
import wildmeshing as wm

def tetrahedralize_bad_mesh(in_file_name, out_file_body="", epsilon=0.001):
    wm.tetrahedralize(in_file_name, out_file_body, epsilon)

def test():
    print("succ")