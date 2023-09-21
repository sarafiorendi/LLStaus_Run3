# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
from nis import match
from unittest import result
import pepper
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb
import coffea
import os

from coffea.nanoevents import NanoAODSchema

class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigSTau
    
    def zero_handler(func):
        def _function(self, data, *args, **kwargs):
            if len(data) > 0: return func(self, data, *args, **kwargs)
            else: return ak.Array([])
        return _function
    
    def __init__(self, config, eventdir):
        # Initialize the class, maybe overwrite some config variables and
        # load additional files if needed
        # Can set and modify configuration here as well
        config["histogram_format"] = "root"
        # Need to call parent init to make histograms and such ready
        super().__init__(config, eventdir)

        if "pileup_reweighting" not in config:
            logger.warning("No pileup reweigthing specified")

        if "jet_fake_rate" not in config and len(config["jet_fake_rate"]) == 0:
            self.predict_jet_fakes = False
            logger.warning("No jet fake rate specified")
        else:
            self.predict_jet_fakes = config["predict_yield"]

    def process_selection(self, selector, dsname, is_mc, filler):
        
        # Triggers
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname,
            is_mc,
            self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))

        if is_mc and "pileup_reweighting" in self.config:
            selector.add_cut("Pileup reweighting", partial(
                self.do_pileup_reweighting, dsname))
            
        # HEM 15/16 failure (2018)
        if self.config["year"] == "2018":
            selector.add_cut("HEM_veto", partial(self.HEM_veto, is_mc=is_mc))

        # MET cut
        selector.add_cut("MET", self.MET_cut)
        
        # 2 muons
        selector.set_column("Muon", self.select_muons)
        selector.set_column("Electron", self.select_electrons)
        selector.add_cut("muon_veto", self.muon_veto)
        selector.add_cut("elec_veto", self.elec_veto)
        
        selector.set_column("Jet_select", self.jet_selection)
        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        
        selector.add_cut("has_two_valid_jets", self.has_two_valid_jets)
        selector.set_column("sum_jj", self.sum_jj)
        
        selector.set_column("mt2_j1_j2_MET", self.get_mt2)
        selector.set_multiple_columns(self.missing_energy)
        selector.set_multiple_columns(self.mt_jets)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        
        selector.add_cut("b_tagged_jet_cut", self.b_tagged_jet_cut)
        
        # selector.add_cut("skim_jets", self.skim_jets)
        
        # Tagger part for calculating scale factors
        # Scale factors should be calculated -
        # before cuts on the number of the jets
        selector.set_multiple_columns(self.set_njets_pass)
        if self.config["predict_yield"]:
            selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))
        
        # Selection of the jet is performed only for two leading jets:

        selector.set_column("Jet_pass", self.jet_passed)
        selector.set_column("Jet_tau", partial(self.jet_tau, is_mc=is_mc))
        
        # # 1. Leading jet is selected by pt
        # selector.add_cut("point", self.b_tagged_jet_cut) # dummy cut to make sure correct Jet_select is used
        # selector.set_column("Jet_select", self.jet_lead)
        # selector.set_multiple_columns(self.set_njets_pass)
        # if self.config["predict_yield"]:
        #     selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))
        # selector.add_cut("redefine_jets_lead2pt", self.b_tagged_jet_cut)
        
        # # 2. Leading jet is selected by tagging score
        # selector.add_cut("point2", self.b_tagged_jet_cut) # dummy cut to make sure correct Jet_select is used
        # selector.set_column("Jet_select", self.jet_lead_score)
        # selector.set_multiple_columns(self.set_njets_pass)
        # if self.config["predict_yield"]:
        #     selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))
        # selector.add_cut("redefine_jets_lead2prob", self.b_tagged_jet_cut)
    
    @zero_handler
    def skim_jets(self, data):
        jets = data["Jet_select"]
        file_name = os.path.basename(data.metadata["filename"])
        dataset_name = data.metadata["dataset"]
        prefix = f"evn_{data.metadata['entrystart']}_{data.metadata['entrystop']}_"
        path = self.config["skim_path"]
        path = f"{path}/{dataset_name}"
        print(path)
        print(f"{path}/"+prefix+file_name+".parquet")
        if not os.path.exists(path):
            os.makedirs(path)
        ak.to_parquet(jets, f"{path}/"+prefix+file_name+".parquet")
        return np.ones(len(data))
    
    @zero_handler
    def HEM_veto(self, data, is_mc):
        weight = np.ones(len(data), dtype=np.float32)
        jets = data["Jet"]
        elctron = data["Electron"]
        electron_in15or16_hem = ( (elctron.pt > 20) & (elctron.eta > -3.0) & (elctron.eta < -1.3) & (elctron.phi > -1.57) & (elctron.phi < -0.87) )
        jet_in15or16_hem = ( (jets.pt > 20) & (jets.eta > -3.2) & (jets.eta < -1.3) & (jets.phi > -1.77) & (jets.phi < -0.67) )
        in_hem = (ak.any(electron_in15or16_hem, axis=-1) | ak.any(jet_in15or16_hem, axis=-1))
        if is_mc:
            weight[in_hem] = (1-0.66)
        else:
            weight[in_hem] = 0.0
        return weight   
    
    @zero_handler
    def jet_lead(self, data):
        # select only first two pt leading jets
        return data["Jet_select"][:,:2]

    @zero_handler
    def jet_passed(self, data):
        jets = data["Jet_select"]
        jets = jets[jets.disTauTag_score1 >= 0.9972]
        return jets
    
    @zero_handler
    def jet_tau(self, data, is_mc):
        jets = data["Jet_pass"]
        if not is_mc: return jets
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_h, dRlist = jets.nearest(tau_vis, return_metric=True, threshold=0.3)
        tau_jet = jets[~ak.is_none(matches_h, axis=1)]
        return tau_jet
    
    @zero_handler
    def jet_lead_score(self, data):
        # select only first two biggest score jets
        jets = data["Jet_select"]
        sort_idx = ak.argsort(jets.disTauTag_score1, axis=-1, ascending=False)
        jets = jets[sort_idx]
        return jets[:,:2]
    
    @zero_handler    
    def get_mt2(self, data):
        jet1 = data["Jet_select"][:,0]
        jet2 = data["Jet_select"][:,1]
        met = data["MET"]
        return mt2.mt2(
            jet1.mass, jet1.px, jet1.py,
            jet2.mass, jet2.px, jet2.py,
            met.px, met.py,
            0, 0
        )
        
    @zero_handler
    def missing_energy(self, data):
        jets = data["Jet_select"]
        HT_valid = ak.sum(jets.pt, axis=-1)
    
        px = ak.sum(jets.px, axis=-1)
        py = ak.sum(jets.py, axis=-1)
        HT_miss_valid = np.sqrt(px*px + py*py)

        jets = data["Jet"]
        HT = ak.sum(jets.pt, axis=-1)

        px = ak.sum(jets.px, axis=-1)
        py = ak.sum(jets.py, axis=-1)
        HT_miss= np.sqrt(px*px + py*py)
        
        return {
            "HT_valid" : HT_valid, "HT_miss_valid" : HT_miss_valid,
            "HT" : HT, "HT_miss" : HT_miss
        }
    
    def delta_phi(self, phi1_ak, phi2_ak):
        phi1 = np.array(phi1_ak)
        phi2 = np.array(phi2_ak)
        assert phi1.shape == phi2.shape
        d = phi1 - phi2
        indx_pos = d>np.pi
        d[indx_pos] -= np.pi*2
        indx_neg = d<=-np.pi
        d[indx_neg] += np.pi*2
        return d
    
    @zero_handler
    def mt_jets(self, data):
        jet1 = data["Jet_select"][:,0]
        jet2 = data["Jet_select"][:,1]
        
        MET = data["MET"]
        one_min_cs = 1.0 - np.cos(self.delta_phi(jet1.phi, MET.phi))
        prod = 2*jet1.pt*MET.pt
        mt_j1 = np.sqrt( prod * one_min_cs)
        
        one_min_cs = 1.0 - np.cos(self.delta_phi(jet2.phi, MET.phi))
        prod = 2*jet2.pt*MET.pt
        mt_j2 = np.sqrt( prod * one_min_cs) 
    
        return {
            "mt_jet1" : mt_j1,
            "mt_jet2" : mt_j2,
            "mt_sum" : mt_j1 + mt_j2
        }
        
    @zero_handler
    def dphi_jet1_jet2(self, data):
        return self.delta_phi(data["Jet_select"][:,0].phi,
                              data["Jet_select"][:,1].phi)
    
    @zero_handler
    def MET_cut(self, data):
        return data["MET"].pt > self.config["MET_cut"]
    
    @zero_handler
    def select_muons(self, data):
        muons = data["Muon"]
        is_good = (
              (muons.pt > self.config["muon_pt_min"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            & (muons[self.config["muon_ID"]] == 1)
            & (muons.pfIsoId >= self.config["muon_pfIsoId"])
            )
        return muons[is_good]
    
    @zero_handler
    def select_electrons(self, data):
        ele = data["Electron"]
        ele_low_eta_iso = eval(self.config["ele_low_eta_iso"])
        ele_high_eta_iso = eval(self.config["ele_high_eta_iso"])
        isolation_cut = ( ele_low_eta_iso | ele_high_eta_iso )
        is_good = (
            isolation_cut
            & (ele.pt > self.config["ele_pt_min"])
            & (ele.eta < self.config["ele_eta_max"])
            & (ele.eta > self.config["ele_eta_min"])
            & (ele[self.config["eleVeto"]] == 1)
            & (ele[self.config["eleID"]] == 1)
            )
        return ele[is_good]
    
    @zero_handler
    def muon_veto(self, data):
        return ak.num(data["Muon"]) == 0
    
    @zero_handler
    def elec_veto(self, data):
        return ak.num(data["Electron"]) == 0

    @zero_handler
    def jet_selection(self, data):
        jets = data["Jet"]
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            & (jets.jetId >= self.config["jet_jetId"] )
            )]
        return jets
    
    @zero_handler
    def has_two_valid_jets(self, data):
        return ak.num(data["Jet_select"]) >= 2
    
    @zero_handler
    def pfcand_valid(self, data):
        pfCands = data["PFCandidate"]
        is_good = (
            (pfCands.pt > self.config["pfcand_pt"])
            & (pfCands.eta < self.config["pfcand_eta_max"])
            & (pfCands.eta > self.config["pfcand_eta_min"])
            & (pfCands[self.config["track"]])
        )
        pfCands_selected = pfCands[is_good]
        sort_idx = ak.argsort(pfCands_selected.pt, axis=-1, ascending=False)
        return pfCands_selected[sort_idx]
    
    @zero_handler
    def match_jet_to_pfcand(self, data, jet_name = None, pf_name = None, dR = 0.4):
        '''
        This function match all particle-flow candidates to every jet in the [jet_name] collection
        (function is alternative to match_nearest, but works for multydim case and expected to be much slower)
        '''
        jets = data[jet_name]
        pfcands = data[pf_name]
        # here return_combinations=True is needed to return _pfcands_unzipped
        # which broadcasted in the way to duplicate pdcands list per every jet in event
        (dr, (_, _pfcands_unzipped)) = jets.metric_table(pfcands, metric=coffea.nanoevents.methods.vector.LorentzVector.delta_r, return_combinations=True)
        pfcands_matched = _pfcands_unzipped[(dr < dR)]
        return pfcands_matched
    
    @zero_handler
    def get_matched_pfCands(self, data, match_object, dR=0.4):
        pfCands = self.match_jet_to_pfcand(data, jet_name=match_object, pf_name="PfCands", dR=dR)
        pfCands_lead = ak.firsts(pfCands, axis=-1)
        pfCands_lead["dxysig"] = pfCands_lead.dxy / pfCands_lead.dxyError
        pfCands_lead["Lrel"] = np.sqrt(pfCands_lead.dxy**2 + pfCands_lead.dz**2)
        pfCands_lead["dxy_weight"] = ak.mean(pfCands.dxy, weight=pfCands.pt, axis=-1)
        pfCands_lead["dxysig_weight"] = ak.mean(pfCands.dxy / pfCands.dxyError, weight=pfCands.pt, axis=-1)
        return pfCands_lead
    
    @zero_handler
    def set_jet_dxy(self, data):
        jets = data["Jet_select"]
        # Mask jets with dxy nan (no selected pfcands matching)
        bad_jets = ak.is_none(data["Jet_lead_pfcand"].dxy, axis=-1)
        jets = ak.mask(jets, ~bad_jets) # mask bad jets to keep coorect shape
        jets["dz"] = np.abs(data["Jet_lead_pfcand"].dz)
        jets["dxy"] = np.abs(data["Jet_lead_pfcand"].dxy)
        jets["dxy_weight"] = np.abs(data["Jet_lead_pfcand"].dxy_weight)
        jets["dxysig"] = np.abs(data["Jet_lead_pfcand"].dxysig)
        jets["dxysig_weight"] = np.abs(data["Jet_lead_pfcand"].dxysig_weight)
        jets = jets[~bad_jets] # remove bad jets
        return jets
    
    @zero_handler
    def sum_jj(self, data):
        return data['Jet_select'][:,0].add(data['Jet_select'][:,1])
    
    @zero_handler
    def b_tagged_jet_cut(self, data):
        # To leading jets are excluded! 
        jet_not_signal = data["Jet_select"][:,2:]
        # Jet_btagDeepFlavB satisfies the Medium (>0.2783) WP:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        b_tagged_idx = (jet_not_signal.btagDeepFlavB > 0.2783)
        return ak.num(jet_not_signal[b_tagged_idx]) == 0
    
    @zero_handler
    def set_njets_pass(self, data):
        jets_score = data["Jet_select"].disTauTag_score1
        n_pass = []
        for score in self.config["score_pass"]:
            jets_pass = jets_score[(jets_score>score)]
            passed = ak.num(jets_pass, axis=1)
            n_pass.append( passed )
        n_pass = ak.from_regular(
            np.stack(n_pass, axis=1), axis=-1)
        return {
            "n_pass" : n_pass,
            "n_pass_score_bin" : ak.local_index(n_pass, axis=1)
        }
        
    @zero_handler
    def predict_yield(self, data, weight=None):
        jets = data["Jet_select"]
        
        # from bin 0 to bin 1
        weights_bin0to1 = []
        for score in self.config["score_pass"][1:]: # skip first bin because it is just 1
            events_0tag = (ak.num(jets[(jets.disTauTag_score1 > score)]) == 0) # events with 0 tag
            jets_notag = (jets.disTauTag_score1 < score) # to have a per jet mask
            jets_counted = jets[events_0tag * jets_notag] # to select only jets in events with 0 tag
            # fake_sf =  self.config["jet_fake_rate"](jet_dxy=jets_counted.dxy, jet_pt=jets_counted.pt, jet_score=score)
            fake_sf =  self.config["jet_fake_rate"](jet_pt=jets_counted.pt, jet_score=score)
            weight_sfs = ak.sum(fake_sf, axis=1)
            weights_bin0to1.append(weight_sfs)
        yield_bin0to1 = ak.from_regular(np.stack(weights_bin0to1, axis=1), axis=-1)
        # print(yield_bin0to1)

        # from bin 1 to bin 2
        weights_bin1to2 = []
        for score in self.config["score_pass"][1:]: # skip first bin because it is just 1
            events_1tag = (ak.num(jets[(jets.disTauTag_score1 > score)]) == 1) # events with 1 tag
            jets_notag = (jets.disTauTag_score1 < score) # to have a per jet mask and not to count the tagged jet
            jets_counted = jets[events_1tag * jets_notag]  # to select only jets in events with 1 tag
            # fake_sf =  self.config["jet_fake_rate"](jet_dxy=jets_counted.dxy, jet_pt=jets_counted.pt, jet_score=score)
            fake_sf =  self.config["jet_fake_rate"](jet_pt=jets_counted.pt, jet_score=score)
            weight_sfs = ak.sum(fake_sf, axis=1)
            weights_bin1to2.append(weight_sfs)
        yield_bin1to2 = ak.from_regular(np.stack(weights_bin1to2, axis=1), axis=-1)
        # print(yield_bin1to2)
        
        # from bin 0 to bin 2
        weights_bin0to2 = []
        for score in self.config["score_pass"][1:]: # skip first bin because it is just 1
            events_0tag = (ak.num(jets[(jets.disTauTag_score1 > score)]) == 0) # events with 0 tag
            jets_notag = (jets.disTauTag_score1 < score) # to have a per jet mask
            jets_counted = jets[events_0tag * jets_notag] # to select only jets in events with 1 tag
            # fake_sf =  self.config["jet_fake_rate"](jet_dxy=jets_counted.dxy, jet_pt=jets_counted.pt, jet_score=score)
            fake_sf =  self.config["jet_fake_rate"](jet_pt=jets_counted.pt, jet_score=score)
            combinations = ak.combinations(fake_sf, 2, axis=1) # to have all possible combinations of 2 jets
            combinations_unzipped = ak.unzip(combinations)
            products = combinations_unzipped[0] * combinations_unzipped[1]
            weight_sfs = ak.sum(products, axis=1)
            weights_bin0to2.append(weight_sfs)
        yield_bin0to2 = ak.from_regular(np.stack(weights_bin0to2, axis=1), axis=-1)
        # print(yield_bin0to2)
        
        # now we need to each predicted yield assign cooresponding score bin
        score_bin = ak.local_index(yield_bin0to1, axis=1) + 1 # +1 because we skip first bin
        # print(score_bin)
        
        return {"yield_bin0to1" : weight*yield_bin0to1,
                "yield_bin1to2" : weight*yield_bin1to2,
                "yield_bin0to2" : weight*yield_bin0to2,
                "score_bin"     : score_bin}
