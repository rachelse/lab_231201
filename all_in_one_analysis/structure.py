import os
from constant import Color as clr
import interface

import nglview as nv
import py3Dmol

class Structure:
    def __init__(self, path = "", file = "", name = ""):
        self.path = path
        self.file = file
        self.name = name

        self.structure = None
        self.chains = None
        self.interface = None
        self.score = None

    def get_structure(self):
        if not self.file.endswith(".pdb"):
            print("File is not a PDB file")
            return None
        path = os.path.join(self.path, self.file) 
        structure = interface.get_structure(name = self.name, path = path)
        chain_name = list(structure[0].child_dict.keys())
        chains =  {}
        for c in chain_name:
            chains[c] = len(structure[0][c])

        self.structure = structure
        self.chains = chains

    def get_interface(self, threshold = 8):
        if self.structure == None:
            print("Predicted structure not loaded")
        elif len(self.chains) != 2:
            print("Predicted structure does not have 2 chains")
            return None

        self.interface = interface.get_interface(self.structure)
    
    def add_score(self, ipTM, iPAE, pdockq):
        self.score = {'ipTM' : ipTM, 'iPAE' : iPAE, 'pdockq' : pdockq}
        return

    def show_3d(self,interface, label, color_scheme, module, monomer):

        if self.structure == None:
            print("Predicted structure not loaded. Try get_structure() first.")
            return None
        
        if len(self.chains) == 1 or monomer:
            if not monomer:
                monomer = list(self.chains.keys())[0]
            str3d = self.set_3dmol_monomer(label, monomer)

        elif len(self.chains) == 2:    
            if module == "nglview":
                str3d = self.set_nglview(interface)
            elif module == "py3dmol":
                str3d = self.set_3dmol(interface, label, color_scheme)

        return str3d

    def set_3dmol(self, interface, label, color_scheme):
        view = py3Dmol.view(height=400, width=700)
        pdbpath = os.path.join(self.path, self.file)
        view.addModel(open(pdbpath,'r').read() , "pdb")
        view.zoomTo()
        
        c1 = list(self.chains.keys())[0]
        c2 = list(self.chains.keys())[1]

        if not interface:
            opacity = 1
        else: opacity = 0.5

        """ Set color scheme """
        if color_scheme == "dpi":
            view.setStyle({'chain': c1}, {'cartoon': {'color': clr.IDR, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c2}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "asyn":
            p1, p2, p3 = 60, 95, 140
            view.setStyle({'chain': c1, 'resi': [*range(1,p1+1)]}, 
                          {'cartoon': {'color': clr.ASYN1, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c1, 'resi': [*range(p1+1,p2+1)]}, 
                          {'cartoon': {'color': clr.ASYN2, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c1, 'resi': [*range(p2+1,p3+1)]}, 
                          {'cartoon': {'color': clr.ASYN3, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c2}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "asyn_p1":
            view.setStyle({'chain': c1}, 
                          {'cartoon': {'color': clr.ASYN1, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c2}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "asyn_p2":
            view.setStyle({'chain': c1}, 
                          {'cartoon': {'color': clr.ASYN2, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c2}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "asyn_p3":
            view.setStyle({'chain': c1}, 
                          {'cartoon': {'color': clr.ASYN3, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain': c2}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "chain":
            view.setStyle({'cartoon': {'colorscheme':'chain', 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "ss":
            view.setStyle({'cartoon': {'colorscheme':'ssJmol', 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "plddt":
            plddt_scheme = {'prop': 'b', 'gradient': 'roygb', 'min': 0, 'max': 100}
            view.setStyle({'cartoon': {'colorscheme': plddt_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        elif color_scheme == "plddt_dpi":
            plddt_scheme = {'prop': 'b', 'gradient': 'roygb', 'min': 0, 'max': 100}
            custom_scheme = {'prop': 'b', 'min': 0, 'max': 100, 'colors':['orange', 'yellow', 'light_blue', 'blue'], 'gradient': 'linear'}
            verylow_scheme = {'prop': 'b', 'min': 0, 'max': 50, 'colors':['orange'], 'gradient': 'linear'}
            low_scheme = {'prop': 'b', 'min': 50, 'max': 70, 'colors':['yellow'], 'gradient': 'linear'}
            confident_scheme = {'prop': 'b', 'min': 70, 'max': 90, 'colors':['light_blue'], 'gradient': 'linear'}
            high_scheme = {'prop': 'b', 'min': 90, 'max': 100, 'colors':['blue'], 'gradient': 'linear'}
   
            # view.setStyle({'chain':c1}, {'cartoon': {'colorscheme': custom_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            view.setStyle({'chain':c1}, {'cartoon': {'colorscheme': plddt_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}})

            # view.setStyle({'chain':c1}, {'cartoon': {'colorfunc': [lambda x: 'orange' if x < 50 else 'yellow' if x < 70 else 'light_blue' if x < 90 else 'blue'], 'tubes': False, 'arrows':True, 'opacity':opacity}})
            # view.addStyle({'chain':c1}, {'cartoon': {'colorscheme': verylow_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}}, 
            #               {'not': {'chain':c1, 'b': [*range(50,100)]}})
            # view.addStyle({'chain':c1}, {'cartoon': {'colorscheme': low_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}},
            #               {'not': {'chain':c1, 'b': [*range(0,50)]}})
            # view.addStyle({'chain':c1}, {'cartoon': {'colorscheme': confident_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            # view.addStyle({'chain':c1}, {'cartoon': {'colorscheme': high_scheme, 'tubes': False, 'arrows':True, 'opacity':opacity}},
            #               {'not': {'chain':c1, 'b': [*range(0,70)]}})

            view.setStyle({'chain': c2}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
        # elif color_scheme == "residue":
        #     # TODO
        #     # blrd = {'prop': 'resi', 'gradient': 'bluered', 'min': 1, 'max': 100}
        #     view.setStyle({'chain': c1},
        #                   {'cartoon': {'color':'spectrum', 'tubes': False, 'arrows':True, 'opacity':opacity}})
        #     view.setStyle({'chain': c2},
        #                  {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':opacity}})
            
        
        # Specify interface by changing opacity
        if interface:
            if color_scheme == "dpi":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]}, 
                              {'cartoon': {'color': clr.IDR, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain':c2, 'resi': self.interface[c2]}, 
                            {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}}) 
            elif color_scheme == "asyn":
                view.addStyle({'chain': c1, 'resi': [i for i in self.interface[c1] if i <= p1]}, 
                              {'cartoon': {'color': clr.ASYN1, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain': c1, 'resi': [i for i in self.interface[c1] if i > p1 and i <= p2]},
                                {'cartoon': {'color': clr.ASYN2, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain': c1, 'resi': [i for i in self.interface[c1] if i > p2 and i <= p3]},
                                {'cartoon': {'color': clr.ASYN3, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain':c2, 'resi': self.interface[c2]}, 
                            {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}}) 
            elif color_scheme == "asyn_p1":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]},
                                {'cartoon': {'color': clr.ASYN1, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain':c2, 'resi': self.interface[c2]},
                            {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}})
            elif color_scheme == "asyn_p2":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]},
                                {'cartoon': {'color': clr.ASYN2, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain':c2, 'resi': self.interface[c2]},
                            {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}})
            elif color_scheme == "asyn_p3":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]},
                                {'cartoon': {'color': clr.ASYN3, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain':c2, 'resi': self.interface[c2]},
                            {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}})
            elif color_scheme == "chain":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]}, {'cartoon': {'colorscheme':'chain', 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain': c2, 'resi': self.interface[c2]}, {'cartoon': {'colorscheme':'chain', 'tubes': False, 'arrows':True, 'opacity':1}})
            elif color_scheme == "ss":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]}, {'cartoon': {'colorscheme':'ssJmol', 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain': c2, 'resi': self.interface[c2]}, {'cartoon': {'colorscheme':'ssJmol', 'tubes': False, 'arrows':True, 'opacity':1}})
            elif color_scheme == "plddt":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]}, {'cartoon': {'colorscheme':plddt_scheme, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain': c2, 'resi': self.interface[c2]}, {'cartoon': {'colorscheme':plddt_scheme, 'tubes': False, 'arrows':True, 'opacity':1}})
            elif color_scheme == "plddt_dpi":
                view.addStyle({'chain': c1, 'resi': self.interface[c1]}, {'cartoon': {'colorscheme':custom_scheme, 'tubes': False, 'arrows':True, 'opacity':1}})
                view.addStyle({'chain': c2, 'resi': self.interface[c2]}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}})
            # elif color_scheme == 'residue':
            #     # TODO
            #     # blrd = {'prop': 'resi', 'gradient': 'bluered', 'min': 1, 'max': 100}
            #     view.addStyle({'chain': c1, 'resi': self.interface[c1]}, {'cartoon': { 'tubes': False, 'arrows':True, 'opacity':1}})
            #     view.addStyle({'chain': c2, 'resi': self.interface[c2]}, {'cartoon': {'color': clr.PARTNER, 'tubes': False, 'arrows':True, 'opacity':1}})

        
        """ Add labels """
        if label:
            cfg = {'backgroundColor': 'white', 'showBackground': False, 'fontColor': 'black', 'fontSize': 12, 'fontFace': 'Helvetica',
                     'fixed':True, 'useScreen':True, 'screenOffset':{'x':15, 'y':0}, 'inFront':True}
            view.addLabel(self.name, cfg)
            cfg['screenOffset']['y'] -= 20
            view.addLabel(f"Chain {c1} {len(self.interface[c1])}/{self.chains[c1]}", cfg)
            cfg['screenOffset']['y'] -= 20
            view.addLabel(f"Chain {c2} {len(self.interface[c2])}/{self.chains[c2]}", cfg)

            if self.score != None:
                cfg['screenOffset']['y'] -= 20
                view.addLabel(f"ipTM{self.score['ipTM']} iPAEScore{self.score['iPAE']} pdockq{self.score['pdockq']}", cfg)
        view.show()
        # return view

    def set_nglview(self, interface):
        if interface and self.interface == None:
            print("Interface not defined")
            return None
        
        if interface:
            interface_scheme = "" #TODO
        ss_scheme = "sstruc"
        
        view = nv.show_biopython(self.structure)

        return view
    
    def set_3dmol_monomer(self, label, chain = "A"):
        view = py3Dmol.view(height=400, width=400)
        path = os.path.join(self.path, self.file)
        view.addModel(open(path,'r').read() , "pdb")
        view.zoomTo()
        return view
    
    def get_color_scheme(self, interface = False):
        #TODO
        return
    
    def save_image():
        #TODO
        return