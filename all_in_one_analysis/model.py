import structure
import information
from score import pdockq, iPAE
import matplotlib.pyplot as plt
import matplotlib.image as img
import json, os

class Model:
    pwd, name, path = None, None, None
    models = None
    
    def __init__(self, pwd, name, path):
        self.pwd = pwd
        self.name = name
        self.path = path

        return

    def get_model(self, rank=1, item = "all"):
        rank = self._get_rank_item(rank)
        self._initialize_model_dic(ranks=rank)
        item = self._get_score_item(item)
        self.get_structure(rank)
        if rank != [""]:
            self.get_score(rank, item)
        return
    
    def get_score(self, rank, item):
        for r in rank:
            if self.models[f'rank{r}']['score'] == None:
                self.models[f'rank{r}']['score'] = {'plddt': None, 'pae': None, 'max_pae': None,'ptm': None, 'iptm': None,
                                                    'pdockq': None, 'iPAE': None}

            
            with open(os.path.join(self.pwd, self.path["score"][f"rank{r}"])) as f:
                score = json.load(f)
                for (k,v) in score.items():
                    if k in item:
                        self.models[f"rank{r}"]['score'][k] = v

            # Other scores not provided by AlphaFold2
            if 'pdockq' in item:
                self.models[f'rank{r}']['score']['pdockq'] = pdockq(os.path.join(self.pwd + self.path["pdb"][f"rank{r}"]))
            if 'iPAE' in item:
                # Structure must be loaded
                if self.models[f'rank{r}']['structure'] == None:
                    self.get_structure(rank=r)
                strc = self.models[f'rank{r}']['structure']
                strc_chain = sorted(list(strc.chains.keys()))
                self.models[f'rank{r}']['score']['iPAE'] = iPAE(strc.interface['pairs'], 
                                                                score['pae'], lenA = strc.chains[strc_chain[0]])
        print(f"Loaded scores of rank{r}")
        print("Find scores by calling e.g. model.models['rank1']['score']['plddt']")

    def get_structure(self, rank):
        for r in rank:
            if self.models[f'rank{r}']['structure'] != None:
                print(f"Structure for rank {r} already loaded")
                continue
            else:
                struc_r = structure.Structure(self.pwd, self.path["pdb"][f"rank{r}"], self.name)
                struc_r.get_structure()
                struc_r.get_interface()
                self.models[f'rank{r}']['structure'] = struc_r
        return 
    
    def get_info(self, uniprotID = ""):
        uniprotID = self._split_name(uniprotID)
        print("Getting information of",uniprotID)
        info = {}
        info['gene'] = self.get_gene(uniprotID)
        info['description'] = self.get_description(uniprotID)
        info['GO'] = self.get_goterms(uniprotID)

        comment = self.get_comment(uniprotID)
        for (k,v) in comment.items():
            info[k] = v
        return info
    
    def get_gene(self, uniprotID = ""):
        uniprotID = self._split_name(uniprotID)
        print("Getting gene name of",uniprotID)
        gene = information.get_gene(uniprotID)
        return gene
    
    def get_goterms(self, uniprotID = "", go = ""):
        uniprotID = self._split_name(uniprotID)
        print("Getting GO terms of",uniprotID)
        goterms = information.get_goterms(uniprotID)
        if go == "id":
            return goterms["id"]
        elif go == "bp":
            return goterms["biological process"]
        elif go == "cc":
            return goterms["cellular component"]
        elif go == "mf":
            return goterms["molecular function"]

        return goterms
    
    def get_description(self, uniprotID = ""):
        uniprotID = self._split_name(uniprotID)
        print("Getting description of",uniprotID)
        return information.get_description(uniprotID)
    
    def get_comment(self, uniprotID = ""):
        uniprotID = self._split_name(uniprotID)
        print("Getting comment of",uniprotID)
        return information.get_comment(uniprotID)
    
    def show_3d(self, rank, interface = False, module = "py3dmol", monomer = False, label = True, color_scheme = "dpi"):

        strc = self.models[f'rank{rank}']['structure']
        try:
            strc.add_score(
                        self.models[f'rank{rank}']['score']['iptm'],
                        self.models[f'rank{rank}']['score']['iPAE'], 
                        self.models[f'rank{rank}']['score']['pdockq'])
        except:
            pass
        strc_3d = strc.show_3d(interface, label, color_scheme, module, monomer)
        return strc_3d
    
    def show_plddt(self):
        return self.show_png("plddt")

    def show_pae(self):
        return self.show_png("pae")
    
    def show_coverage(self):
        return self.show_png("coverage")
        
    def show_png(self, item):
        # item : pae, plddt, coverage
        item = item + ".png"
        pngimg = img.imread(os.path.join(self.pwd, self.path[item]))
        plt.imshow(pngimg)
        plt.axis('off')
        return plt.show()
    
    def _get_rank_item(self, rank):
        if rank == "all":
            rank = self.path["score"].keys()
        elif rank == "":
            rank = [""]
        else:
            if type(rank) == str:
                rank = [int(rank)]
            elif type(rank) == list:
                rank = [int(i) for i in rank]
            elif type(rank) == int:
                rank = [rank]
        return rank

    def _get_score_item(self, item = "all"):
        if item == "all":
            item = ["plddt", "pae", "ptm","max_pae", "iptm",'pdockq', 'iPAE']
        else:
            if type(item) == str:
                item = [item]
            elif type(item) == list:
                item = item
        return item
    
    def _initialize_model_dic(self, ranks):
        init_dic = {'score': None, 
                    'structure': None}
        self.models = {}
        for r in ranks:
            self.models[f'rank{r}'] = init_dic.copy()

    def _split_name(self, query ):
        if query in ["partner", ""]:
            return self.name.split("_")[1]
        elif query in ["target", "idr"]:
            return self.name.split("_")[0]
        return query