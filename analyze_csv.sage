from sage.features.latte import Latte_count, Latte_integrate
from matplotlib.pyplot import figure
from sage.graphs.connectivity import connected_components_subgraphs
from numpy import concatenate as sumarr
import csv
import os 
if "flow_meanders.sage.py" in os.listdir():
    if "fm.py" in os.listdir():
        os.remove("fm.py")
    os.rename("flow_meanders.sage.py", "fm.py")
import fm

def params_global():
    params_global = ["No.", "top_comp", "bottom_comp", "num_pieces", "cycles", "paths", "points", "index",
        "fm_dimension", "fm_volume", "fm_nvolume", "fm_latpts", "fm_ehrhart", "fm_fvector", 
        "sfm_dimension", "sfm_volume", "sfm_nvolume", "sfm_latpts", "sfm_ehrhart", "sfm_fvector"]
    return params_global

"""
FILL OUT REST OF DOCSTRING
vecs: Can be a single vector or a list of vectors. 
      If it is a single vector, we truncate the right side 
      to make it the appropriate size for every n.
"""
def makecsv_allcomps(low, high, vecs, progbar=1000):
    params = params_global()
    counter = 1
    total_pairs = sum([binomial(2^(i-1), 2) for i in range(low, high+1)])
    if isinstance(vecs, tuple):
        vecs = [tuple(vecs[:i]) for i in range(low, high+1)]

    with open("All_Composition_Pairs_with_Weights_" + str(low) + "-" + str(high) + ".csv", 'w') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=params, delimiter=";")
        writer.writeheader()
        for n in range(low, high+1):
            vec = vecs[n-low]
            comps = list(Compositions(n))
            for i in range(len(comps)):
                for j in range(i, len(comps)):
                    data = fm.data_dictionary(list(comps[i]), list(comps[j]), vec)
                    data["No."] = counter 
                    counter += 1
                    writer.writerow(data)
                    if progbar and counter % progbar == 0:
                        print(counter, "/", total_pairs, "done")

def makecsv_ncomps(low, high, vecs, progbar=1000):
    params = params_global()
    counter = 1
    total_pairs = sum([2^(i-1) for i in range(low, high+1)])
    if isinstance(vecs, tuple):
        vecs = [tuple(vecs[:i]) for i in range(low, high+1)]

    with open("All_Compositions_Over_(n)_with_Weights_" + str(low) + "-" + str(high) + ".csv", 'w') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=params, delimiter=";")
        writer.writeheader()
        for n in range(low, high+1):
            vec = vecs[n-low]
            comps = list(Compositions(n))
            for i in range(len(comps)):
                data = fm.data_dictionary(list(comps[i]), [n], vec)
                data["No."] = counter 
                counter += 1
                writer.writerow(data)
                if progbar and counter % progbar == 0:
                    print(counter, "/", total_pairs, "done")

def makecsv_onescomps(low, high, vecs, progbar=1000):
    params = params_global()
    counter = 1
    total_pairs = sum([2^(i-1) for i in range(low, high+1)])
    if isinstance(vecs, tuple):
        vecs = [tuple(vecs[:i]) for i in range(low, high+1)]

    with open("All_Compositions_Over_(1^n)_with_Weights_" + str(low) + "-" + str(high) + ".csv", 'w') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=params, delimiter=";")
        writer.writeheader()
        for n in range(low, high+1):
            vec = vecs[n-low]
            comps = list(Compositions(n))
            for i in range(len(comps)):
                data = fm.data_dictionary(list(comps[i]), [1]*n, vec)
                data["No."] = counter 
                counter += 1
                writer.writerow(data)
                if progbar and counter % progbar == 0:
                    print(counter, "/", total_pairs, "done")

def get_matches(csvfpath, param1, param2):
    filename = "_".join(["matches", param1, param2, ".txt"])
    out = open(filename, "a")
    with open(csvfpath) as csvf:
        reader = csv.DictReader(csvf, delimiter=";")
        for row in reader:
            if row[param1] == row[param2]:
                print(row["top_comp"], "//", row["bottom_comp"], "\n", param1, row[param1], "\n", param2, row[param2], "\n", file=out)
    out.close()

makecsv_ncomps(2,12,(1,0,0,0,0,0,0,0,0,0,0,0))
makecsv_onescomps(2,12,(1,0,0,0,0,0,0,0,0,0,0,0))
# compare = ["fm_volume", "fm_nvolume", "fm_latpts", "sfm_volume", "sfm_nvolume", "sfm_latpts"]
# for param in compare:
#     get_matches("All_Composition_Pairs_with_Weight_2-8.csv", "index", param)