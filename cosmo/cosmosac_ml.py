import numpy as np
import os
import pandas as pd
from cosmosac import CosmoSac
import util as u

# keras version = 2.6.0 nightly
from keras.models import load_model

modelpath = u.here + "\\model\\"


class CosmoSacGcm(CosmoSac):
    file = pd.ExcelFile(modelpath + "\\Regressed MG.xlsx")
    vol = pd.read_excel(file, sheet_name="Area Volume", index_col=0)["Volume"]
    sig = pd.read_excel(file, sheet_name="Sigma Profile", index_col=0)
    nhb = pd.read_excel(file, sheet_name="NHB Sigma Profile", index_col=0)
    oh = pd.read_excel(file, sheet_name="OH Sigma Profile", index_col=0)
    ot = pd.read_excel(file, sheet_name="OT Sigma Profile", index_col=0)

    def add_comp(self, smiles, name=None):
        first = u.mol_frag(smiles, "MG First", 3)
        second = u.mol_frag(smiles, "MG Second", duplicate=True)
        mg = pd.concat([first, second])

        V = 0
        sigma_profiles = np.zeros((3, 51))

        for frag_num in mg.index:
            V += mg.loc[frag_num] * self.vol.loc[frag_num]
            sigma_profiles[0] += mg.loc[frag_num] * self.nhb.loc[frag_num]
            sigma_profiles[1] += mg.loc[frag_num] * self.oh.loc[frag_num]
            sigma_profiles[2] += mg.loc[frag_num] * self.ot.loc[frag_num]

        sigma_profiles = np.where(sigma_profiles < 0, 0, sigma_profiles)
        A = np.sum(sigma_profiles)

        self.A.append(A)
        self.V.append(V)
        self.psigA.append(sigma_profiles)
        self.name.append(name)

    def get_totsig(self, smiles):
        first = u.mol_frag(smiles, "MG First", 3)
        second = u.mol_frag(smiles, "MG Second", duplicate=True)
        mg = pd.concat([first, second])

        psigA = 0

        for frag_num in mg.index:
            psigA += mg.loc[frag_num] * self.sig.loc[frag_num]
        psigA = np.where(psigA < 0, 0, psigA)

        return psigA

    def gam(self):
        ln_gam_comb = self.ln_gam_comb()
        ln_gam_res = self.ln_gam_res()

        ln_gam = ln_gam_comb + ln_gam_res
        gam = np.exp(ln_gam)
        return gam


class CosmoSacMgcn(CosmoSac):
    sig = load_model(modelpath + "mgcn_sig.h5")
    nhb = load_model(modelpath + "mgcn_nhb.h5")
    oh = load_model(modelpath + "mgcn_oh.h5")
    ot = load_model(modelpath + "mgcn_ot.h5")
    vol = load_model(modelpath + "mgcn_vol.h5")

    def add_comp(self, smiles, name=None):
        nfm, efm = u.mgcn_input(smiles)
        nfm = np.reshape(nfm, [1, *np.shape(nfm)])
        efm = np.reshape(efm, [1, *np.shape(efm)])

        # area = np.sum(self.sig.predict(pred_input))
        volume = 759 * self.vol.predict([nfm, efm])[0][0]
        sigma_profiles = np.zeros((3, 51))
        sigma_profiles[0] = 202 * self.nhb.predict([nfm, efm])[0]
        sigma_profiles[1] = 7 * self.oh.predict([nfm, efm])[0]
        sigma_profiles[2] = 16 * self.ot.predict([nfm, efm])[0]
        sigma_profiles = np.where(sigma_profiles < 0, 0, sigma_profiles)
        area = np.sum(sigma_profiles)

        self.A.append(area)
        self.V.append(volume)
        self.psigA.append(sigma_profiles)
        self.name.append(name)

    def get_totsig(self, smiles):
        nfm, efm = u.mgcn_input(smiles)
        nfm = np.reshape(nfm, [1, *np.shape(nfm)])
        efm = np.reshape(efm, [1, *np.shape(efm)])

        psigA = 202 * self.sig.predict([nfm, efm])
        psigA = np.where(psigA < 0, 0, psigA)

        return psigA

    def gam(self):
        ln_gam_comb = self.ln_gam_comb()
        ln_gam_res = self.ln_gam_res()

        ln_gam = ln_gam_comb + ln_gam_res
        gam = np.exp(ln_gam)
        return gam


class CosmoSacGcgcn(CosmoSac):
    sig = load_model(modelpath + "gcgcn_sig.h5")
    nhb = load_model(modelpath + "gcgcn_nhb.h5")
    oh = load_model(modelpath + "gcgcn_oh.h5")
    ot = load_model(modelpath + "gcgcn_ot.h5")
    vol = load_model(modelpath + "gcgcn_vol.h5")

    def add_comp(self, smiles, name=None):
        nfm, efm = u.gcgcn_input(smiles)
        nfm = np.reshape(nfm, [1, *np.shape(nfm)])
        efm = np.reshape(efm, [1, *np.shape(efm)])

        # area = np.sum(self.sig.predict(mat_input))
        volume = 562 * self.vol.predict([nfm, efm])[0][0]
        sigma_profiles = np.zeros((3, 51))
        sigma_profiles[0] = 145 * self.nhb.predict([nfm, efm])[0]
        sigma_profiles[1] = 7 * self.oh.predict([nfm, efm])[0]
        sigma_profiles[2] = 16 * self.ot.predict([nfm, efm])[0]
        sigma_profiles = np.where(sigma_profiles < 0, 0, sigma_profiles)
        area = np.sum(sigma_profiles)

        self.A.append(area)
        self.V.append(volume)
        self.psigA.append(sigma_profiles)
        self.name.append(name)

    def get_totsig(self, smiles):
        nfm, efm = u.gcgcn_input(smiles)
        nfm = np.reshape(nfm, [1, *np.shape(nfm)])
        efm = np.reshape(efm, [1, *np.shape(efm)])

        psigA = 145 * self.sig.predict([nfm, efm])[0]
        psigA = np.where(psigA < 0, 0, psigA)

        return psigA

    def gam(self):
        ln_gam_comb = self.ln_gam_comb()
        ln_gam_res = self.ln_gam_res()

        ln_gam = ln_gam_comb + ln_gam_res
        gam = np.exp(ln_gam)
        return gam
