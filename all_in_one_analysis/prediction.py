import os
import matplotlib.pyplot as plt
import matplotlib.image as img

import loader
import model

class Prediction:
    pwd = None
    model_paths = None

    def __init__(self, pwd, item = "all", model = None):
        if not os.path.exists(pwd):
            print("Path does not exist")
            return None

        if type(model) == str:
            model = [model]

        item = self._get_item(item)
            
        self.pwd = pwd
        self.model_paths = self.extract_files(pwd, item, model)
        return

    def extract_files(self, pwd, item, model):
        files = sorted(os.listdir(pwd))
        initial = {i:None for i in item}

        model_paths = {}

        for f in files:
            loaded = loader._parse_model(f)
            if loaded:
                name, i = loaded
                if model and (name not in model):
                    continue
                elif i not in item:
                    continue
                elif item == "done":
                    continue

                if name not in model_paths.keys():
                    model_paths[name] = initial.copy()

                if i == "score" or i == "pdb":
                    if model_paths[name][i] == None:
                        model_paths[name][i] = {}
                    try:
                        r = int(loader._parse_rank(f))
                    except:
                        r = ""
                    model_paths[name][i][f"rank{r}"] = f
                else: model_paths[name][i] = f

        return model_paths
    
    def show_png(self, model, item ):
        # item : pae, plddt, coverage
        item = item + ".png"
        pngimg = img.imread(os.path.join(self.pwd, self.model_paths[model][item]))
        plt.imshow(pngimg)
        plt.axis('off')

    def load_model(self, name, item = "all"):
        item = self._get_item(item)
        info = {}
        for i in item:
            info[i] = self.model_paths[name][i]
        m = model.Model(self.pwd, name, info)
        return m
    
    def _get_item(self, item="all"):
        if item == "all":
            item = ["msa", "pae.png", "plddt.png", "coverage.png", "pae", "score", "pdb"]
        else:
            if type(item) == str:
                item = [item]
            elif type(item) == list:
                item = item
        return item


