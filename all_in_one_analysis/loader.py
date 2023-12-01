import os, re

def get_path(work):
    work = work.lower()
    work_path = {
        'benchmark' : "/home/seamustard52/dpi/benchmark/prediction_dimer/",
        'benchmark_old': "/home/seamustard52/dpi_old/benchmark/prediction_dimer_v2/",
        'benchmark_exp' : "/home/seamustard52/dpi_old/benchmark/pdb",
        'control' : "/home/seamustard52/dpi_old/benchmark/prediction_control/",
        'asyn' : "/home/seamustard52/dpi_old/asyn_screen/prediction/",
        'asyn_p1': "/home/seamustard52/dpi_old/asyn_screen/prediction_p1/",
        'asyn_p2': "/home/seamustard52/dpi_old/asyn_screen/prediction_p2/",
        'asyn_p3': "/home/seamustard52/dpi_old/asyn_screen/prediction_p3/",
        'brca' : "/home/seamustard52/dpi_old/brca2_screen/prediction/",
        'granzyme' : "/home/seamustard52/granzymeE/prediction/"
    }
    return work_path[work]

def _parse_model(fname):
    if fname.endswith(".done.txt") or fname == "config.json":
        return False
    elif fname.endswith(".a3m"):
        model = fname.split(".")[0]
        item = "msa"
    elif fname.endswith(".png"):
        model = '_'.join(fname.split(".")[0].split("_")[:-1]) 
        item = fname.split("_")[-1].lower() # pae.png or plddt.png or coverage.png 
    elif fname.endswith(".pdb"):
        model = re.findall(r"(.+)_unrelaxed_rank_", fname)
        if len(model) == 0:
            model = re.findall(r"(.+).pdb", fname)
        model = model[0]
        item = "pdb"
    elif fname.endswith(".json"):
        if fname.endswith("_predicted_aligned_error_v1.json"):
            model = re.findall(r"(.+)_predicted_aligned_error", fname)[0]
            item = "pae"
        else:
            model = '_'.join(re.findall(r"(.+)_scores", fname)[0].split("_"))
            item = "score"
    else:
        return False
    return model, item

def _parse_rank(fname):
    if fname.endswith('.pdb'):
        rank = re.findall(r"rank_(\d+)_", fname)[0]
    elif fname.endswith('.json') and "scores" in fname:
        rank = re.findall(r"rank_(\d+)_", fname)[0]
    else:
        print("rank not found")
        return None
    return int(rank)

def get_item(item, version="v3"):
    item_expr = {
        'msa': {'v2': '.a3m', 'v3':'.a3m'},
        'pae_plot' : {'v2': "_PAE.png", 'v3':"_pae.png"},
        'plddt': {'v2': '_plddt.png', 'v3':"_plddt.png"},
        'coverage': {'v2': '_coverage.png', 'v3':'_coverage.png'},
        'pae' :{'v2': '_predicted_aligned_error', 'v3':'_predicted_aligned_error'},
        'score' : {'v2': '_scores', 'v3':'_scores'},
        'pdb' : {'v2': '.pdb', 'v3':'.pdb'}
    }
    return item_expr[item][version]

def get_all_files(work, item = "all", version="v3", path = ""):
    # work: benchmark, control, asyn, brca
    # item: msa, pae_plot, plddt, coverage, pae, score, pdb
    # version: v2, v3
    if len(path) > 0:
        work_dir = path
    else:
        work_dir = get_path(work)
    if item == "all":
        files = os.listdir(work_dir)
    else:
        files = [f for f in os.listdir(work_dir) if get_item(item, version) in f]
    return files

def get_model_files(model, work, item="all", version="v3", path = ""):
    # model: 
    # work: benchmark, control, asyn, brca
    # item: msa, pae_plot, plddt, coverage, pae, score, pdb
    # version: v2, v3
    item = item.lower()
    if work == "asyn":
        model = f"P37840_{model}"
    elif work == "brca":
        model = f"P51587_{model}"
    print(f"Finding {item} files of {model}")

    if len(path) >0:
        work_dir = path
    else:
        work_dir = get_path(work) 

    files = sorted(os.listdir(work_dir))
    model_files = []
    for f in files:
        if item == "all":
            if f.startswith(model):
                print(f)
                model_files.append(f)
        else:
            pattern = get_item(item, version)
            if f.startswith(model) and pattern in f:
                print(f)
                model_files.append(f)

    return model_files