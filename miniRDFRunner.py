#!/usr/bin/env python3

from python.mkShapesRDFMini import shapeRDFMaker
import ROOT
ROOT.EnableImplicitMT()
import numpy as np
import uproot
import time

cuts = ['res_sig_ele_DNN_geq_0p5',
 'res_sig_mu_DNN_geq_0p5',
 'res_sig_ele_DNN_leq_0p5',
 'res_sig_mu_DNN_leq_0p5',
 'boost_sig_ele_DNN_geq_0p5',
 'boost_sig_mu_DNN_geq_0p5',
 'boost_sig_ele_DNN_leq_0p5',
 'boost_sig_mu_DNN_leq_0p5',
 'res_wjetcr_ele',
 'res_wjetcr_mu',
 'boost_wjetcr_ele',
 'boost_wjetcr_mu',
 'res_topcr_ele',
 'res_topcr_mu',
 'boost_topcr_ele',
 'boost_topcr_mu']

res_cuts = [i for i in cuts if i.startswith("res_")]
boost_cuts = [i for i in cuts if i.startswith("boost_")]
sig_cuts = [i for i in cuts if "sig" in i]

variables = {}

variables['events']  = {   'name': 'events',      
                        'range' : (1,0,2),  
                        'xaxis' : 'events', 
                        'fold' : 3
                        }

########################


variables['Mww'] = {   'name': 'Mww',
                        'range' : ([ 200., 400., 600., 800., 1000., 1200., 1500., 2000., 3000.],), #variable range  
                        'xaxis' : 'M_{WV} [GeV]',
                        'fold' : 3,
                        'blind': [200. ,3000.]
                        }

variables['Mww_resolved_6'] = {   'name': 'Mww',
                        'range' : ([ 200., 400., 600., 800., 1000., 1200., 2000],), #variable range  
                        'xaxis' : 'M_{WV} [GeV]',
                        'fold' : 3,
                        'blind': [200. ,2000.]
                        }


#####################
#Fit variables

variables['fit_bins_res'] ={  'name' : 'fit_bin_res',
                            'range' : (21,1,22),
                            'xaxis' : 'Wjets resolved bin', 
                            'fold' : 0,
                            'cuts': res_cuts
}   

variables['fit_bins_boost'] ={  'name' : 'w_lep_pt',
                            'range' : ([0,50,100,150,200,300,400,600],),
                            'xaxis' : 'W leptonic Pt', 
                            'fold' : 3,
                            'cuts': boost_cuts
}   



pdf_var =  ['pdf_weight_1718V{}Var'.format(i) for i in range(103)]

full = {k: "Resolved" for k in cuts if k.startswith("res_")}
boost = {k: "Boosted" for k in cuts if k.startswith("boost_")}

full.update(boost)

sm = shapeRDFMaker("../PlotsConfigurations/Configurations/VBSjjlnu/Full2018v7_EFT/rootFile_fit_2018_EFT_DNNcut_geq0p5_sm_2into4_trees/plots_fit_2018_EFT_DNNcut_geq0p5_sm_2into4_trees_ALL_sm.root")
sm.defineSample("sm") 
sm.defineVariables(variables)
sm.defineStructure(
    full
)

subs = {
        "sm": "weight_sm",
        "lin_cHW": "weight_lin_cHW",
        "quad_cHW": "weight_quad_cHW",
        "lin_cW": "weight_lin_cW",
        "quad_cW": "weight_quad_cW",
        "lin_cHWB": "weight_lin_cHWB",
        "quad_cHWB": "weight_quad_cHWB",
        "lin_cHbox": "weight_lin_cHbox",
        "quad_cHbox": "weight_quad_cHbox",
        "lin_cHQ1": "weight_lin_cHQ1",
        "quad_cHQ1": "weight_quad_cHQ1",
        "lin_cHj1": "weight_lin_cHj1",
        "quad_cHj1": "weight_quad_cHj1",
        "lin_cHl1": "weight_lin_cHl1",
        "quad_cHl1": "weight_quad_cHl1",
        "lin_cHB": "weight_lin_cHB",
        "quad_cHB": "weight_quad_cHB",
        "linear_mixed_cHWB_cHW": "weight_linear_mixed_cHWB_cHW",
        "linear_mixed_cHB_cHj1": "weight_linear_mixed_cHB_cHj1",
        "linear_mixed_cHWB_cHbox": "weight_linear_mixed_cHWB_cHbox",
        "linear_mixed_cW_cHWB": "weight_linear_mixed_cW_cHWB",
        "linear_mixed_cW_cHB": "weight_linear_mixed_cW_cHB",
        "linear_mixed_cHWB_cHB": "weight_linear_mixed_cHWB_cHB",
        "linear_mixed_cHbox_cHl1": "weight_linear_mixed_cHbox_cHl1",
        "linear_mixed_cW_cHW": "weight_linear_mixed_cW_cHW",
        "linear_mixed_cHbox_cHB": "weight_linear_mixed_cHbox_cHB",
        "linear_mixed_cHWB_cHl1": "weight_linear_mixed_cHWB_cHl1",
        "linear_mixed_cHl1_cHj1": "weight_linear_mixed_cHl1_cHj1",
        "linear_mixed_cHW_cHl1": "weight_linear_mixed_cHW_cHl1",
        "linear_mixed_cW_cHl1": "weight_linear_mixed_cW_cHl1",
        "linear_mixed_cHW_cHB": "weight_linear_mixed_cHW_cHB",
        "linear_mixed_cW_cHj1": "weight_linear_mixed_cW_cHj1",
        "linear_mixed_cHW_cHj1": "weight_linear_mixed_cHW_cHj1",
        "linear_mixed_cW_cHbox": "weight_linear_mixed_cW_cHbox",
        "linear_mixed_cHbox_cHW": "weight_linear_mixed_cHbox_cHW",
        "linear_mixed_cHB_cHQ1": "weight_linear_mixed_cHB_cHQ1",
        "linear_mixed_cW_cHQ1": "weight_linear_mixed_cW_cHQ1",
        "linear_mixed_cHl1_cHB": "weight_linear_mixed_cHl1_cHB",
        "linear_mixed_cHWB_cHQ1": "weight_linear_mixed_cHWB_cHQ1",
        "linear_mixed_cHWB_cHj1": "weight_linear_mixed_cHWB_cHj1",
        "linear_mixed_cHbox_cHQ1": "weight_linear_mixed_cHbox_cHQ1",
        "linear_mixed_cHbox_cHj1": "weight_linear_mixed_cHbox_cHj1",
        "linear_mixed_cHQ1_cHj1": "weight_linear_mixed_cHQ1_cHj1",
        "linear_mixed_cHW_cHQ1": "weight_linear_mixed_cHW_cHQ1",
        "linear_mixed_cHl1_cHQ1": "weight_linear_mixed_cHl1_cHQ1"
        
    }

# subs = {
#         "sm": "weight_sm",
#         "lin_cHW": "weight_lin_cHW",
#         "quad_cHW": "weight_quad_cHW",
#     }

sm.defineSubSamples(
    subs
)
sm.excludeNuisSubsamples(
    {k: pdf_var for k in subs.keys()}
)

sm.defineEnvelopes(
    {
           "pdf_weight_1718": {
               "samples": ["sm"],
               "variation": pdf_var
           }
    }
)

nuisAlias = {"sm": {"CMS_PS_FSRDown": "CMS_PS_FSR_VBSDown", "CMS_PS_FSRUp": "CMS_PS_FSR_VBSUp", "CMS_PS_ISRDown": "CMS_PS_ISR_VBSDown", "CMS_PS_ISRUp": "CMS_PS_ISR_VBSUp"}}
nuisAlias.update({k: {"CMS_PS_ISRDown": "CMS_PS_ISR_EFTDown", "CMS_PS_ISRUp": "CMS_PS_ISR_EFTUp", "CMS_PS_FSRDown": "CMS_PS_FSR_EFTDown", "CMS_PS_FSRUp": "CMS_PS_FSR_EFTUp"} for k in subs.keys()})

sm.defineNuisanceAlias(nuisAlias)

start_time = time.time()
print("--- %s seconds ---" % (time.time() - start_time))

sm.buildRDFs()

sm.buildHistograms()


sm.convertHistoToNumpy()

sm.computeEnvelopes()

print("--- %s seconds ---" % (time.time() - start_time))

sm.writeToFile("prova2.root")
