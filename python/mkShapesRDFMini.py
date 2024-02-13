import ROOT
ROOT.EnableImplicitMT()
import numpy as np
import uproot

class shapeRDFMaker:
    df_dict = {}
    histo_numpy_dict = {}
    nuisance_reweight = {}
    nuisance_tree = {}
    histo_dict = {}
    filepath = ""
    subsamples = {}
    sample = ""
    file = ""
    exclude_nuisances_subs = {}
    envelopes = {}
    histo_list = []
    structure = {}
    variables = {}
    nuisanceAlias = {}    

    def __init__(self, file):
        self.file = file
    
    @staticmethod
    def __getHistoProto(name_, range_):
        if len(range_) == 1:
            edges = np.array(range_[0], dtype=np.double)
            return (name_, name_, len(edges)-1, edges)
        elif len(range_) == 3:
            return (name_, name_, *range_)
        else:
            raise ValueError("Histo protos should have len 1 ([...], ) for variable bin lengths and length 3 (bins, min, max) for equal size bins")
        
    def __mergeStructure(self):
        
        # self.structure has keys the cuts and values the tree names
        # while self.variables has the variable definition need to merge the two 
        
        self.structure__ = {}
        for region, tree in self.structure.items():
            self.structure__[region] = {
                "tree": tree,
                "variables": []
            }
            for aliasvariable, vdefinition in self.variables.items():
                if "cuts" in vdefinition and region not in vdefinition["cuts"]: continue
                vdefinition["alias"] = aliasvariable
                self.structure__[region]["variables"].append(vdefinition)
        
        self.structure = self.structure__
        
    def defineNuisanceAlias(self, d):
        self.nuisanceAlias = d
        
    def defineVariables(self, d):
        self.variables = d
        if self.structure: self.__mergeStructure()
        
    def defineStructure(self, d):
        self.structure = d
        if self.variables: self.__mergeStructure()
    
    def defineSample(self, sample):
        self.sample = sample
    
    def defineSubSamples(self, d):
        self.subsamples = d
    
    def defineEnvelopes(self, d):
        self.envelopes = d
        
    def excludeNuisSubsamples(self, excludeNSub):
        self.exclude_nuisances_subs = excludeNSub
        
    def writeToFile(self, file_, target="histo_dict"):
        
        target = getattr(self, target)
        file = uproot.recreate(file_)
        for r in target.keys():
            for v in target[r].keys():
                for h in target[r][v].keys():
                    try:
                        file[r + "/" + v + "/" + h] = target[r][v][h]
                    except TypeError:
                        file[r + "/" + v + "/" + h] = target[r][v][h].GetValue()
 
    def computeEnvelopes(self):
        
        if not self.histo_numpy_dict:
            print("---> Will first convert TH1D to Numpy")
            if not self.histo_dict:
                print("Need to run shapes first")
                return 
            self.convertHistoToNumpy()
        
        for envelopeName, envelopeIngredients in self.envelopes.items():
            for r in self.histo_numpy_dict.keys():
                for v in self.histo_numpy_dict[r].keys():
                    histo_names = self.histo_numpy_dict[r][v].keys()
                    for sample in envelopeIngredients["samples"]:
                        upName = "histo_" + sample + "_" + envelopeName + "Up"
                        downName = "histo_" + sample + "_" + envelopeName + "Down"
                        if upName in histo_names or downName in histo_names: continue
                        if not all(i in histo_names for i in ["histo_" + sample + "_" + j for j in envelopeIngredients["variation"]]): 
                            print("Not all ingredients present, will skip")
                            continue
                            
                        # just get bin edges, they are all the same
                        be__ = self.histo_numpy_dict[r][v]["histo_" + sample + "_" + envelopeIngredients["variation"][0]][1]
                            
                        concatHist = tuple(np.resize(self.histo_numpy_dict[r][v][i][0], (1,self.histo_numpy_dict[r][v][i][0].shape[0])) for i in ["histo_" + sample + "_" + j for j in envelopeIngredients["variation"]])
                        UpVar, DownVar = np.amax(concatHist, axis=0), np.amin(concatHist, axis=0)
                        
                        #reshapeback
                        UpVar = np.resize(UpVar, (UpVar.shape[1],))
                        DownVar = np.resize(DownVar, (DownVar.shape[1],))
                        
                        # also save these as TH1D
                        hu = ROOT.TH1D(r + "_" + v + "_" + upName, upName, len(be__)-1, be__)
                        hd = ROOT.TH1D(r + "_" + v + "_" + downName, downName, len(be__)-1, be__)
                        
                        for i in range(len(be__)-1):
                            hu.SetBinContent(i+1, UpVar[i])
                            hd.SetBinContent(i+1, DownVar[i])
                            
                        self.histo_dict[r][v][upName] = hu
                        self.histo_dict[r][v][downName] = hd
                        
                        self.histo_numpy_dict[r][v][upName] = (UpVar, be__)
                        self.histo_numpy_dict[r][v][downName] = (DownVar, be__)
        
    
    @staticmethod
    def TH1ToNumpy(histo):
        histo = histo.GetValue()
        vals = np.array(histo)
        ranges = np.array([histo.GetBinLowEdge(i+1) for i in range(histo.GetNbinsX()+1)])
        # resize the vector to exclude underflow and overflow
        tv = vals[1:-1]
        assert len(ranges)-1 == len(tv), print("Attenzione")
        
        return (tv, ranges)
        
    def convertHistoToNumpy(self):
        print("Starting")
        for region in self.histo_dict.keys():
            self.histo_numpy_dict[region] = {}
            for variable in self.histo_dict[region].keys():
                self.histo_numpy_dict[region][variable] = {}
                for hname in self.histo_dict[region][variable].keys():
                    histo = self.TH1ToNumpy(self.histo_dict[region][variable][hname])
                    self.histo_numpy_dict[region][variable][hname] = histo
                    
        print("Done")
    
    def buildRDFs(self):
        
        f = ROOT.TFile(self.file)
        regions = self.structure.keys()
        for r in regions:
            self.df_dict[r] = {}
            d = f.Get(r + "/" + self.structure[r]["tree"])
            trees = [i.GetName() for i in d.GetListOfKeys()]

            nominal = "tree_" + self.sample
            # this trees will have 2 branches one named "weight" and the other named as "var"
            varied_trees = {k: j for k,j in [[i,i.split(nominal + "_")[1]] for i in trees if i != nominal]}

            # 
            nom_t = ROOT.TChain(r + "/" + self.structure[r]["tree"] + "/" + nominal)
            nom_t.Add(self.file)
            nom_t.SetDirectory(0)

            # Store reweight
            rew_branches = [j.split("reweight_")[1] for j in [i.GetName() for i in nom_t.GetListOfBranches() if i.GetName().startswith("reweight_")]]

            #
            for variation, alias in varied_trees.items():
                vt = ROOT.TChain(r + "/" + self.structure[r]["tree"] + "/" + variation)
                vt.Add(self.file)
                vt.SetDirectory(0)
                nom_t.AddFriend(vt, alias)

            df = ROOT.RDataFrame(nom_t)
            # df = ROOT.RDF.Experimental.Distributed.Spark.RDataFrame(nom_t)
            # define events variable
            df = df.Define("events", "1")

            self.df_dict[r]["df"] = df
            self.df_dict[r]["reweight"] = rew_branches
            self.df_dict[r]["trees"] = [j for _,j in varied_trees.items()]
            self.df_dict[r]["chain"] = nom_t
    
    def runTheHisto(self):
        print("Running graphs")
        # print(histo_list)
        ROOT.RDF.RunGraphs(self.histo_list)
        
    def foldHistos(self):
        for region in self.histo_dict.keys():
            for variable in self.histo_dict[region].keys():
                if "fold" in self.variables[variable]:
                    doFold = self.variables[variable]["fold"]
                else:
                    continue
                    
                for hname in self.histo_dict[region][variable].keys():
                    if doFold == 2 or doFold == 3 :
                       # Setting underflow to first bin
                       uf = self.histo_dict[region][variable][hname].GetBinContent(0)
                       fb = self.histo_dict[region][variable][hname].GetBinContent(1)
                       
                       uf_err = self.histo_dict[region][variable][hname].GetBinError(0)
                       fb_err = self.histo_dict[region][variable][hname].GetBinError(1)
                       
                       self.histo_dict[region][variable][hname].SetBinContent(1, uf + fb)
                       self.histo_dict[region][variable][hname].SetBinError(1, uf_err + fb_err)
                       self.histo_dict[region][variable][hname].SetBinContent(0, 0)
                       self.histo_dict[region][variable][hname].SetBinError(0, 0)
                       
                       # Setting Overflow to last bin
                       nbins = self.histo_dict[region][variable][hname].GetNbinsX()
                       
                       of = self.histo_dict[region][variable][hname].GetBinContent(nbins + 1)
                       lb = self.histo_dict[region][variable][hname].GetBinContent(nbins)
                       
                       of_err = self.histo_dict[region][variable][hname].GetBinError(nbins + 1)
                       lb_err = self.histo_dict[region][variable][hname].GetBinError(nbins)
                       
                       self.histo_dict[region][variable][hname].SetBinContent(nbins, of + lb)
                       self.histo_dict[region][variable][hname].SetBinError(nbins, of_err + lb_err)
                       self.histo_dict[region][variable][hname].SetBinContent(nbins + 1, 0)
                       self.histo_dict[region][variable][hname].SetBinError(nbins + 1, 0)                
        print("Done")
        
                    
    def buildHistograms(self):
        
        for region in self.structure.keys():
            print(region)
            self.histo_dict[region] = {}
            for variable in self.structure[region]["variables"]:
                # print(variable)
                 
                alias = variable["alias"]
                tree_name = variable["name"]
                range_ = variable["range"]
                fold = variable["fold"]
                
                self.histo_dict[region][alias] = {}
                
                # nominal histogtram
                nominal = "histo_" + self.sample
                histo_proto = self.__getHistoProto(nominal, range_)
                self.histo_dict[region][alias][nominal] = self.df_dict[region]["df"].Histo1D(histo_proto, tree_name, "weight")
                self.histo_list.append(self.histo_dict[region][alias][nominal])
                
                # print("Done nominal")
                
                
                # work on tree-type variations
                rew_var = self.df_dict[region]["reweight"]
                for rv in rew_var:
                    new_w = "reweight_" + rv
                    nname = rv
                    if self.sample in self.nuisanceAlias and rv in self.nuisanceAlias[self.sample]: nname = self.nuisanceAlias[self.sample][rv]
                    var_name = "histo_" + self.sample + "_" + nname
                    histo_proto = self.__getHistoProto(var_name, range_)
                    self.histo_dict[region][alias][var_name] = self.df_dict[region]["df"].Define("ovweight", "weight*{}".format(new_w)).Histo1D(histo_proto, tree_name, "ovweight")
                    self.histo_list.append(self.histo_dict[region][alias][var_name])
                        
                #print("Done reweight variation")
                # work on suffix type variations
                tree_var = self.df_dict[region]["trees"]
                for tv in tree_var:
                    new_w = tv + ".weight"
                    nname = tv 
                    if self.sample in self.nuisanceAlias and tv in self.nuisanceAlias[self.sample]: nname = self.nuisanceAlias[self.sample][tv]
                    var_name = "histo_" + self.sample + "_" + nname
                    histo_proto = self.__getHistoProto(var_name, range_)
                    if tv + "." + tree_name in self.df_dict[region]["df"].GetColumnNames():
                       self.histo_dict[region][alias][var_name] = self.df_dict[region]["df"].Define("ovweight", new_w).Histo1D(histo_proto, tv + "." + tree_name, "ovweight")
                    else:
                       # print("Column {} not present in df, will use {} instead".format(tv + "." + tree_name, tree_name))
                       self.histo_dict[region][alias][var_name] = self.df_dict[region]["df"].Define("ovweight", new_w).Histo1D(histo_proto, tree_name, "ovweight")

                    self.histo_list.append(self.histo_dict[region][alias][var_name])
                    
                #print("Done tree variation")
                
                # do the same for the subsamples
                for sName, sWeight in self.subsamples.items():
                    
                    # nominal histogtram
                    nominal = "histo_" + sName
                    histo_proto = self.__getHistoProto(nominal, range_)
                    self.histo_dict[region][alias][nominal] = self.df_dict[region]["df"].Define("ovweight", f"weight*{sWeight}").Histo1D(histo_proto, tree_name, "ovweight")
                    self.histo_list.append(self.histo_dict[region][alias][nominal])
                    
                    # work on tree-type variations
                    rew_var = [ i for i in self.df_dict[region]["reweight"] if sName in self.exclude_nuisances_subs and i not in self.exclude_nuisances_subs[sName] ]
                    for rv in rew_var:
                        new_w = f"{sWeight}*reweight_" + rv
                        nname = rv
                        if self.sample in self.nuisanceAlias and rv in self.nuisanceAlias[self.sample]: nname = self.nuisanceAlias[self.sample][rv]
                        var_name = "histo_" + sName + "_" + nname
                        histo_proto = self.__getHistoProto(var_name, range_)
                        self.histo_dict[region][alias][var_name] = self.df_dict[region]["df"].Define("ovweight", "weight*{}".format(new_w)).Histo1D(histo_proto, tree_name, "ovweight")
                        self.histo_list.append(self.histo_dict[region][alias][var_name])
                            
                    #print("Done reweight variation")
                    # work on suffix type variations
                    tree_var = [ i for i in self.df_dict[region]["trees"] if sName in self.exclude_nuisances_subs and i not in self.exclude_nuisances_subs[sName] ]
                    for tv in tree_var:
                        new_w = tv + ".weight*" + tv + "." + f"{sWeight}"
                        nname = tv
                        if self.sample in self.nuisanceAlias and tv in self.nuisanceAlias[self.sample]: nname = self.nuisanceAlias[self.sample][tv]
                        var_name = "histo_" + sName + "_" + nname
                        histo_proto = self.__getHistoProto(var_name, range_)
                        if tv + "." + tree_name in self.df_dict[region]["df"].GetColumnNames():
                           self.histo_dict[region][alias][var_name] = self.df_dict[region]["df"].Define("ovweight", new_w).Histo1D(histo_proto, tv + "." + tree_name, "ovweight")
                        else:
                           # print("Column {} not present in df, will use {} instead".format(tv + "." + tree_name, tree_name))
                           self.histo_dict[region][alias][var_name] = self.df_dict[region]["df"].Define("ovweight", new_w).Histo1D(histo_proto, tree_name, "ovweight")

                        self.histo_list.append(self.histo_dict[region][alias][var_name])
        
        print("Finished histogram booking")
        # run histos in multiprocess
        self.runTheHisto()
        
        print("Finished filling")
        print("Folding")
        
        self.foldHistos()
                    
