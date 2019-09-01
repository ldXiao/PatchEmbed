import igl
import pybcclean as pbc
import polyfempy as pf
import wildmeshing as wm
import meshplot as mp
import os
import numpy as np
from scipy.spatial import KDTree
import yaml
def barycenters(v, f):
    s = np.zeros((len(f),3),dtype='float64')
    for i, r in enumerate(s):
        s[i] = (v[f[i][0]]+v[f[i][1]]+v[f[i][2]])/3
    return s

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
    count =0
    for label, surface in enumerate(yaml_dict["surfaces"]):
        for face_id in surface["face_indices"]:
            face_label_dict[face_id] = label
            count +=1

    with open(out_file_rel_dir, 'w') as f:
        f.write("1"+" "+str(count)+"\n")
        for item in sorted(face_label_dict.keys()):
            f.write("%s\n" % face_label_dict[item])
    return out_file_rel_dir

def separate_surfaces(v, f, fl):
    """
    :param v vertices
    :param f faces
    :param fl face labels
    perturb the surfaces according to the face labels
    :return v_dict, f_dict
    """
    v_dict = {}
    f_dict_temp = {}
    f_dict = {}
    count_dict = {}
    for fidx, lb in enumerate(fl):
        if lb not in f_dict_temp:
            f_dict_temp[lb] = np.zeros_like(f)
            f_dict_temp[lb][0,:]=f[fidx,:].copy()
            count_dict[lb] = 1
        else:
            f_dict_temp[lb][count_dict[lb],:]=f[fidx,:].copy()
            count_dict[lb] += 1
    for lb in f_dict_temp:
        f_dict_temp[lb] = f_dict_temp[lb][0:count_dict[lb],:]
        v_dict[lb], f_dict[lb], _, _= igl.remove_unreferenced(v,f_dict_temp[lb])
    return v_dict, f_dict

def perturb_and_union(v_dict, f_dict, eps):
    """
    :eps float to control the perturbation
    return new_v, new_f, new_fl
    """
    nv = 0 
    nf = 0
    for lb in f_dict:
        nv += len(v_dict[lb])
        nf += len(f_dict[lb])
    new_v = np.zeros((nv,3),dtype=float)
    new_f = np.zeros((nf,3), dtype=int)
    new_fl = np.zeros((nf,1), dtype=int)
    count_v = 0
    count_f =0
    for lb in f_dict:
        f_temp = f_dict[lb].copy()
        v_temp = v_dict[lb].copy()
        new_v[count_v:count_v+len(v_temp)]= v_temp + np.ones_like(v_temp) * np.random.uniform(-eps, eps)
        f_temp += np.ones_like(f_temp) * count_v
        count_v += len(v_temp)
        new_f[count_f:count_f+len(f_temp)]=f_temp
        new_fl[count_f:count_f+len(f_temp)]= np.ones((len(f_temp),1), dtype=int) *lb
        count_f += len(f_temp)
    return new_v, new_f, new_fl

file_bench = "3/00000008_9b3d6a97e8de4aa193b81000_trimesh_000.obj"
file_yml = "3/00000008_9b3d6a97e8de4aa193b81000_features_000.yml"
file_root = os.path.dirname(file_bench)
# wm.tetrahedralize("2/00000006_d4fe04f0f5f84b52bd4f10e4_trimesh_001.obj", file_root+"/"+"bench.mesh", stop_quality=7)
fl_bench_file = parse_feat(file_yml)
fl_bench = igl.read_dmat(fl_bench_file)
v_bench, f_bench = igl.read_triangle_mesh(file_bench)
v_ini, f_ini=igl.read_triangle_mesh("bench.mesh__sf.obj")
prob_mat_ini, fl_ini_temp =pbc.project_face_labels(v_bench,f_bench.astype('int32'), fl_bench.astype('int32'),v_ini,f_ini.astype('int32'))
fl_ini = pbc.refine_labels(v_ini, f_ini.astype('int32'), prob_mat_ini, fl_ini_temp.astype('int32'),1)
# mp.plot(v_ini, f_ini, fl_ini_temp[:,0])
# mp.plot(v_ini, f_ini, fl_ini[:,0])
eps = 0.1
# v_dict, f_dict =separate_surfaces(v_bench, f_bench, fl_bench)
v_dict, f_dict =separate_surfaces(v_ini, f_ini, fl_ini[:,0])
v_bad, f_bad, fl_bad = perturb_and_union(v_dict, f_dict, eps)
# mp.plot(v_bad, f_bad, fl_bad[:,0], shading={"wireframe":True})
eps_file_name = file_root+"/"+str(eps)+"pertb.obj"
igl.write_obj(eps_file_name, v_bad, f_bad)
v_good, f_good = igl.read_triangle_mesh(file_root+"/"+"out.mesh__sf.obj")
print("py0")
print("called")
prob_mat, fl_good_proj = pbc.project_face_labels(v_bad, f_bad.astype('int32'), fl_bad[:,0].astype('int32'), v_good, f_good.astype('int32'))
print("first")
# prob_mat, fl_good_proj = pbc.project_face_labels(v_bad, f_bad.astype('int32'), fl_bad[:,0].astype('int32'), v_good, f_good.astype('int32'))
print("second")
# igl.write_dmat(file_root+"/"+"feat_cut.dmat",fl_good_proj)