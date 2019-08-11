import yaml
import os
import wildmeshing as wm
import numpy as np

def tetrahedralize_bad_mesh(in_obj_file_name, out_file_body="", stop_quality=10):
    parent_dir = os.path.dirname(in_obj_file_name)
    out_file_name = str(stop_quality)+"__sf.obj"
    for file in os.listdir(parent_dir):
        if str(stop_quality)+"_.msh" in file:
            return os.path.join(parent_dir, out_file_name)

    wm.tetrahedralize(in_obj_file_name, os.path.join(parent_dir, str(stop_quality)), stop_quality=stop_quality)
    return os.path.join(parent_dir,out_file_name)

def parse_feat(in_feat_file_name, out_file_body="feat"):
    parent_dir = os.path.dirname(in_feat_file_name)
    out_file_name = out_file_body+".dmat"
    out_file_rel_dir = os.path.join(parent_dir, out_file_name)
        # return the file path if the output file already exists
    with open(in_feat_file_name, 'r') as stream:
        try:
            yaml_dict=yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    face_label_dict = {}
    count =0;
    for label, surface in enumerate(yaml_dict["surfaces"]):
        for face_id in surface["face_indices"]:
            face_label_dict[face_id] = label % 6
            count +=1

    with open(out_file_rel_dir, 'w') as f:
        f.write("1"+" "+str(count)+"\n")
        for item in sorted(face_label_dict.keys()):
            f.write("%s\n" % face_label_dict[item])
    return out_file_rel_dir

def kth_existing_file(grandparent_dir, k):
    count = 0
    parent_dirs = []
    count = -1
    for dir in sorted(os.listdir(grandparent_dir)):
        if "." in dir:
            continue
        if os.listdir(os.path.join(grandparent_dir, dir)):
            parent_dirs.append(dir)
            count += 1
            if count == k:
                break
    par_dir =os.path.join(grandparent_dir, parent_dirs[-1])
    for file in os.listdir(par_dir):
        if ".yml" in file:
            return os.path.join(par_dir, file)
        if (file[-4:]==".obj"):
            if len(file)> 8 and file[-8:]=="__sf.obj":
                continue
            else:
                return os.path.join(par_dir, file)
    else:
        return None


