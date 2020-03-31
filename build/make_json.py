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
    # compute max and min intensities for all of V_full
    vr = V_full.ravel()
    (lvr,) = vr.shape
    v2d = vr.reshape((lvr//3, 3))
    vxy = v2d[:, :2]
    norms_squared = norm(vxy, axis=1)# ** 2
    max_intensity = norms_squared.max()
    min_intensity = norms_squared.min()
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
    w('"min_intensity": %s,\n' % min_intensity)
    w('"max_intensity": %s,\n' % max_intensity)
    inside = False
    for level in range(num_levels):
        if inside:
            w(",\n")
        else:
            w("\n")
        Elevel = level_json(E_full, V_full, level, min_intensity, max_intensity)
        w('"E%s": %s' % (level, Elevel))
        inside = True
    w("}\n")
    out.close()
    return prefix

def level_json(E_full, V_full, level, min_intensity, max_intensity):
    (nl, nr, nc, nlv, three) = V_full.shape
    assert three == 3
    (nl, nr, nc, nlv) = E_full.shape
    E = E_full[:,:,:,level].reshape((nl, nr, nc))
    V = V_full[:,:,:,level,:].reshape((nl, nr, nc, 3))
    #C = colorizer(V)
    xy_colors = xy_colorize(V, min_intensity, max_intensity)
    C = xy_colors["colors"]
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

def ary(*vals):
    return np.array(vals, dtype=np.float)

# tetrahedral colors
vcolors = [
    # vertex, color for tetrahedron
    (ary(-1,-1,-1), ary(0.7,0.7,0.7)),
    (ary(1,1,-1), ary(1,1,0)),
    (ary(1,-1,1), ary(1,0,1)),
    (ary(-1,1,1), ary(0,1,1)),
]

def tetrahedral_colors(velocities, vcolors=vcolors):
    assert velocities.shape[-1] == 3
    vr = velocities.ravel()
    (lvr,)  = vr.shape
    velocities2d = vr.reshape((lvr//3, 3))
    norms = norm(velocities2d, axis=1)
    max_norm = norms.max()
    nvelocities = velocities2d/max_norm
    colors = nvelocities[:]
    for (i, nv) in enumerate(nvelocities):
        nm = norm(nv)
        c = 0
        for (vertex, color) in vcolors:
            c += (vertex.dot(nv)) * color
        c[c < 0] = 0  # null out negatives
        cn = norm(c)
        if cn < 0.001:
            cn = 1.0  # hack
            nm = 0
        colors[i] = (nm/cn) * c
    return colors

# cubic colors
ccolors = [
    # vertex, color for cube
    (ary(1,0,0), ary(1,0,0)),
    (ary(0,1,0), ary(0,1,0)),
    (ary(0,0,1), ary(0,0,1)),
    (ary(-1,0,0), ary(0,.7,.7)),
    (ary(0,-1,0), ary(.7,0,.7)),
    (ary(0,0,-1), ary(.7,.7,0)),
]

verbose = False

def cubic_colors(velocities):
    vcolors = ccolors
    assert velocities.shape[-1] == 3
    vr = velocities.ravel()
    (lvr,)  = vr.shape
    velocities2d = vr.reshape((lvr//3, 3))
    norms = norm(velocities2d, axis=1)
    max_norm = norms.max()
    nvelocities = velocities2d/max_norm
    colors = nvelocities[:]
    for (i, nv) in enumerate(nvelocities):
        nm = norm(nv)
        c = nv * 0.0
        for (vertex, color) in vcolors:
            d = vertex.dot(nv)
            if d > 0:
                c += d * color
        c[c < 0] = 0  # null out negatives
        cn = norm(c)
        if cn < 0.001:
            cn = 1.0  # hack
            nm = 0
        color = (nm/cn) * c
        if verbose:
            print(nv, "colorized", color)
        colors[i] = color
    return colors

#colorizer = tetrahedral_colors
colorizer = cubic_colors

def xy_colorize(velocities, max_intensity, min_intensity):
    max_color = ary(1.0, 0, 0)
    min_color = ary(0.0, 1.0, 1.0)
    #max_intensity = min_intensity = None
    assert velocities.shape[-1] == 3
    vr = velocities.ravel()
    (lvr,)  = vr.shape
    velocities2d = vr.reshape((lvr//3, 3))
    vxy = velocities2d[:, :2]
    norms_squared = norm(vxy, axis=1)# ** 2
    #max_intensity = norms_squared.max()
    #min_intensity = norms_squared.min()
    colors = velocities2d[:]
    diff = max_intensity - min_intensity
    for (i, ns) in enumerate(norms_squared):
        lmbda = (ns - min_intensity) / diff
        colors[i,:] = (1 - lmbda) * min_color + lmbda * max_color
        if ns == min_intensity or ns == max_intensity:
            print (ns, "extreme color", colors[i])
    return {
        "max_color": list(max_color),
        "min_color": list(min_color),
        "max_intensity": max_intensity,
        "min_intensity": min_intensity,
        "colors": colors,
    }

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
    global verbose
    verbose = True
    a = np.arange(2*3*5).reshape((2,3,5)) * 0.5
    j = ravelled_json(a)
    a2 = np.array(json.loads(j)).reshape(a.shape)
    assert(a.tolist() == a2.tolist())
    velocities = np.arange(3*4*5*2*3).reshape((3,4,5,2,3))
    colors = colors_from_velocities(velocities)
    assert colors.shape == velocities.shape
    prefix = jsoniffy("XXX_Fakedata.h5")
    json_path = os.path.join(JSON_FOLDER, prefix+".json")
    f = open(json_path)
    D = json.load(f)
    verbose = False

if __name__ == "__main__":
    #test()
    jsoniffy_all()

    