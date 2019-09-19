import sys
sys.path.append('.')
import utlis as ut 
import wildmeshing as wm
import os

def main():
    granparent_dir = "../data/"
    for num in range(10,15):
        parent_dir = "../data/"+str(num)+"/"
        yaml_file=ut.get_feat_file(parent_dir)
        print("got yml file "+ yaml_file)
        dmat_file =ut.parse_feat(yaml_file)
        print("stored dmat file"+ dmat_file)
        bad_mesh_obj_file = ut.get_bad_mesh(parent_dir)
        if("good.mesh" not in os.listdir(parent_dir)):
            print("tetrhedralizing....")
            if(num != 10):
                wm.tetrahedralize(bad_mesh_obj_file, parent_dir+"good.mesh", stop_quality=7)
            else:
                wm.tetrahedralize(bad_mesh_obj_file, parent_dir+"good.mesh", stop_quality=10)
    return 0

if __name__ == "__main__":
    main()