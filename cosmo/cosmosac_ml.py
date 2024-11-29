import numpy as np
import os
import sys
import pandas as pd

sys.path.append(os.path.dirname(__file__))
from cosmosac import CosmoSac
import util as u

# keras version = 2.6.0 nightly
from keras.models import load_model

modelpath = u.here + "\\model\\"


class CosmoSacGcgcn(CosmoSac):
    sig = load_model(modelpath + "gcgcn_sig.h5")
    nhb = load_model(modelpath + "gcgcn_nhb.h5")
    oh = load_model(modelpath + "gcgcn_oh.h5")
    ot = load_model(modelpath + "gcgcn_ot.h5")
    vol = load_model(modelpath + "gcgcn_vol.h5")

    def add_water(self):
        pass

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
