"""
Construct JSON from raw data in HDF5 files.
"""

import os
import json
import h5py
import numpy as np
from numpy.linalg import norm

RAW_FOLDER = "../raw_data"
JSON_FOLDER = "../docs/json"
CONFIG_PATH = "../docs/config.json"

def jsoniffy_all():
    # remove existing json files
    oldjsons = os.listdir(JSON_FOLDER)
    for oldjson in oldjsons:
        if oldjson.endswith(".json"):
            oldpath = os.path.join(JSON_FOLDER, oldjson)
            print ("removing", oldpath)
            os.unlink(oldpath)
    filenames = os.listdir(RAW_FOLDER)
    prefixes = [jsoniffy(filename) for filename in filenames if filename.endswith(".h5")]
    D = {"prefixes": sorted(prefixes)}
    f = open(CONFIG_PATH, "w")
    json.dump(D, f)
    f.close()
    print("wrote configuration", repr(CONFIG_PATH))

def jsoniffy(filename):
    [prefix, ext] = filename.split(".")
    assert ext == "h5"
    src_path = os.path.join(RAW_FOLDER, filename)
    dst_path = os.path.join(JSON_FOLDER, prefix + ".json")
    print("Converting", repr(src_path), "to json in", repr(dst_path))
    f = h5py.File(src_path, 'r')
    E_full = np.array(f["E_full"], dtype=np.float)
    V_full = np.array(f["V_full"], dtype=np.float)
    (num_layers, num_rows, num_cols, num_levels) = E_full.shape
    assert num_levels >= 2
    (nl, nr, nc, nlv, three) = V_full.shape
    assert (nl, nr, nc, nlv) == (num_layers, num_rows, num_cols, num_levels)
    assert three == 3
    min_value = E_full.min()
    max_value = E_full.max()
    out = open(dst_path, "w")
    w = out.write
    w("{\n")
    w('"num_levels": %s,\n' % num_levels)
    w('"num_layers": %s,\n' % num_layers)
    w('"num_rows": %s,\n' % num_rows)
    w('"num_columns": %s,\n' % num_cols)
    w('"min_value": %s,\n' % min_value)
    w('"max_value": %s,\n' % max_value)
    inside = False
    for level in range(num_levels):
        if inside:
            w(",\n")
        else:
            w("\n")
        Elevel = level_json(E_full, V_full, level)
        w('"E%s": %s' % (level, Elevel))
        inside = True
    w("}\n")
    out.close()
    return prefix

def level_json(E_full, V_full, level):
    (nl, nr, nc, nlv, three) = V_full.shape
    assert three == 3
    (nl, nr, nc, nlv) = E_full.shape
    E = E_full[:,:,:,level].reshape((nl, nr, nc))
    V = V_full[:,:,:,level,:].reshape((nl, nr, nc, 3))
    C = colors_from_velocities(V)
    L = []
    a = L.append
    a("{\n")
    a('"values": %s,' % ravelled_json(E))
    a('"colors": %s' % ravelled_json(C))
    a("}")
    return "".join(L)

def colors_from_velocities(velocities):
    assert velocities.shape[-1] == 3
    vr = velocities.ravel()
    (lvr,)  = vr.shape
    velocities2d = vr.reshape((lvr//3, 3))
    norms = norm(velocities2d, axis=1)
    max_norm = norms.max()
    colors2d = np.abs(velocities2d/max_norm)
    colors = colors2d.reshape(velocities.shape)
    return colors

def float_fmt(x):
    return "%3.2e" % x

def array_content(array):
    "Dump ravelled array content as string."
    s = array.shape
    ls = len(s)
    L = []
    a = L.append
    inside = False
    if len(s) == 1:
        for x in array:
            if inside:
                a(",")
            a(float_fmt(x))
            inside = True
    else:
        assert len(s)>1
        for ar in array:
            if inside:
                a(",\n")
            a(array_content(ar))
            inside = True
    return "".join(L)

def ravelled_json(array):
    return "[\n" + array_content(array) + "\n]"

def test():
    a = np.arange(2*3*5).reshape((2,3,5)) * 0.5
    j = ravelled_json(a)
    a2 = np.array(json.loads(j)).reshape(a.shape)
    assert(a.tolist() == a2.tolist())
    velocities = np.arange(3*4*5*2*3).reshape((3,4,5,2,3))
    colors = colors_from_velocities(velocities)
    assert colors.shape == velocities.shape
    json_path = jsoniffy("Fake.h5")
    f = open(json_path)
    D = json.load(f)

if __name__ == "__main__":
    #test()
    jsoniffy_all()

    