import numpy as np
from scipy.spatial import distance_matrix


class CosmoSac:
    # Some codes are cited from https://doi.org/10.1021/acs.jctc.9b01016


    def __init__(self, version=2013):
        # version and system
        self.ver = version  # COSMO-SAC version
        self.x = []  # liquid mole fraction
        self.T = 0  # system temperature

        # molecular information from COSMO calculation
        self.A = []  # cavity area
        self.V = []  # cavity volume
        self.psigA = []  # sigma profile * area
        self.mtype = []  # molecular dispersive type
        self.e = []  # dispersive parameter /K
        self.name = []  # molecule name

        # global parameters
        self._q0 = 79.53  # area normalization parameter, ang**2
        self._r0 = 66.69  # volume normalization parameter, ang**3
        self._z = 10  # coordination number
        self._sighb = 0.0084  # hydrogen bonding screening charge, e/ang**2
        self._R = 1.987204258e-3  # gas constant, # kcal/K/mol
        self._fdecay = 0.52928 ** (-2)  # unit conversion parameter, 1
        self._sig0 = 0.007  # hydrogen bondable screening charge, e/ang**2
        self._AES = 6525.69  # electrostatic constant A, kcal*ang**4/mol/e**2
        self._BES = 1.4859e8
        # electrostatic constant B, kcal*ang**4*K**2/mol/e**2
        self._aeff, self._n_psig, self._chb, self._cES = self._get_var(version)
        # effective area, number of sigma profiles,
        # hydrogen bonding parameter, electrostatic parameter
        # ang**2, 1, kcal*ang**4/mol/e**2, kcal*ang**4/mol/e**2
        self._reff = np.sqrt(self._aeff / np.pi)  # effective radius, ang
        self._rc = {
            "H": 0.31,
            "He": 0.28,
            "Li": 1.28,
            "Be": 0.96,
            "B": 0.84,
            "C": 0.76,  # sp3 hybridization, sp2: 0.73 sp1: 0.69
            "N": 0.71,
            "O": 0.66,
            "F": 0.57,
            "Ne": 0.58,
            "Na": 1.66,
            "Mg": 1.41,
            "Al": 1.21,
            "Si": 1.11,
            "P": 1.07,
            "S": 1.05,
            "Cl": 1.02,
            "Ar": 1.06,
            "K": 2.03,
            "Ca": 1.76,
            "Sc": 1.70,
            "Ti": 1.60,
            "V": 1.53,
            "Cr": 1.39,
            "Mn": 1.39,  # l.s.; h.s.: 1.61
            "Fe": 1.32,  # l.s.; h.s.: 1.52
            "Co": 1.26,  # l.s.; h.s.: 1.50
            "Ni": 1.24,
            "Cu": 1.32,
            "Zn": 1.22,
            "Ga": 1.22,
            "Ge": 1.20,
            "As": 1.19,
            "Se": 1.20,
            "Br": 1.20,
            "Kr": 1.16,
            "Rb": 2.20,
            "Sr": 1.95,
            "Y": 1.90,
            "Zr": 1.75,
            "Nb": 1.64,
            "Mo": 1.54,
            "Tc": 1.47,
            "Ru": 1.46,
            "Rh": 1.42,
            "Pd": 1.39,
            "Ag": 1.45,
            "Cd": 1.44,
            "In": 1.42,
            "Sn": 1.39,
            "Sb": 1.39,
            "Te": 1.38,
            "I": 1.39,
            "Xe": 1.40,
            "Cs": 2.44,
            "Ba": 2.15,
            "La": 2.07,
            "Ce": 2.04,
            "Pr": 2.03,
            "Nd": 2.01,
            "Pm": 1.99,
            "Sm": 1.98,
            "Eu": 1.98,
            "Gd": 1.96,
            "Tb": 1.94,
            "Dy": 1.92,
            "Ho": 1.92,
            "Er": 1.89,
            "Tm": 1.90,
            "Yb": 1.87,
            "Lu": 1.87,
            "Hf": 1.75,
            "Ta": 1.70,
            "W": 1.62,
            "Re": 1.51,
            "Os": 1.44,
            "Ir": 1.41,
            "Pt": 1.36,
            "Au": 1.36,
            "Hg": 1.32,
            "Tl": 1.45,
            "Pb": 1.46,
            "Bi": 1.48,
            "Po": 1.40,
            "At": 1.50,
            "Rn": 1.50,
            "Fr": 2.60,
            "Ra": 2.21,
            "Ac": 2.15,
            "Th": 2.06,
            "Pa": 2.00,
            "U": 1.96,
            "Np": 1.90,
            "Pu": 1.87,
            "Am": 1.80,
            "Cm": 1.69,
        }
        # atom radius, ang

        # unit change
        self._ang_au = 0.52917721067  # angstrom per 1 atomic unit

    def _get_var(self, ver):
        poss_ver = (2002, 2010, 2013, 2019)

        if ver == 2002:
            aeff = 7.5
            n_psig = 1
            chb = np.zeros((n_psig, n_psig)) + 85580
            cES = lambda T: 8233.36

        elif ver == 2010 or ver == 2013:
            aeff = 7.25
            n_psig = 3  # nhb, OH, OT
            chb = np.zeros((n_psig, n_psig))
            chb[1][1] = 4013.78
            chb[1][2] = chb[2][1] = 3016.43
            chb[2][2] = 932.31
            cES = lambda T: self._AES + self._BES / T / T

        elif ver == 2019:
            aeff = 7.25
            n_psig = 4  # nhb, OH, OT, COOH
            chb = np.zeros((n_psig, n_psig))
            chb[1][1] = 4013.78
            chb[1][2] = chb[2][1] = 3016.43
            chb[1][3] = chb[3][1] = 3020.18
            chb[2][2] = 932.31
            chb[2][3] = chb[3][2] = 1872.84
            chb[3][3] = 2225.67
            cES = lambda T: self._AES + self._BES / T / T

        else:
            raise ValueError(f"The version must be one of {poss_ver}.")

        return aeff, n_psig, chb, cES

    def _is_ms(self, openfile):
        line = openfile.readline()
        if "COSMO Results from DMol3" in line:  # Material Studio 2017
            return True
        elif "text" in line:  # VT database (2006), UD database (2020)
            return False
        else:
            raise ValueError("Initial text in the cosmo file is not valid.")

    def get_cosmo(self, file):
        # reading cosmo file
        file = open(file, "r")
        is_ms = self._is_ms(file)
        flag = "default"

        # cosmo data storage
        atom = []  # atom symbol
        coord = []  # atom coordinate
        seg = []
        # row : segment index
        # col : atom index, xyz positions, area, charge/area

        for line in file:

            if is_ms:

                # flag change
                if "$coordinates xyz [au]" in line:
                    flag = "coordinate"
                    continue
                if "n  atom        position (X, Y, Z) [au]" in line:
                    flag = "segment"
                    continue

                # data extraction
                if "Surface area of cavity" in line:
                    line = line.split()
                    A = float(line[6])  # au**2
                if "Total Volume of cavity" in line:
                    line = line.split()
                    V = float(line[6])  # au**3

                if flag == "coordinate":
                    if "$end" in line:
                        flag = "default"
                    else:
                        line = line.split()
                        atom.append(line[0])
                        coord.append(list(map(float, line[1:4])))  # au
                if flag == "segment":
                    line = line.split()
                    if line:
                        seg.append(
                            [int(line[1]) - 1] + list(map(float, line[2:5] + line[6:8]))
                        )
                        # 0, au, au, au, au**2, e/au**2
            else:

                # flag change
                if "!DATE" in line:
                    flag = "coordinate"
                    continue
                if "n   atom        position (X, Y, Z) [au]" in line:
                    flag = "segment"
                    continue

                # data extraction
                if "Total surface area of cavity" in line:
                    line = line.split()
                    A = float(line[7])  # ang**2
                if "Total volume of cavity" in line:
                    line = line.split()
                    V = float(line[6])  # ang**3

                if flag == "coordinate":
                    if "end" in line:
                        flag = "default"
                    else:
                        line = line.split()
                        atom.append(line[7])
                        coord.append(list(map(float, line[1:4])))  # ang
                if flag == "segment":
                    line = line.split()
                    if line:
                        seg.append(
                            [int(line[1]) - 1] + list(map(float, line[2:5] + line[6:8]))
                        )
                        # 0, au, au, au, ang**2, e/ang**2

        # list to numpy.ndarray
        atom = np.array(atom)
        coord = np.array(coord)
        seg = np.array(seg)

        # unit change
        if is_ms:
            A *= self._ang_au**2
            V *= self._ang_au**3
            coord *= self._ang_au
            seg[:, 1:4] *= self._ang_au
            seg[:, 4] *= self._ang_au**2
            seg[:, 5] /= self._ang_au**2
        else:
            seg[:, 1:4] *= self._ang_au

        return A, V, atom, coord, seg

    def get_bond(self, atom, coord):
        n_atom = len(atom)
        d = distance_matrix(coord, coord)

        matrix = np.zeros((n_atom, n_atom))
        for i in range(n_atom):
            for j in range(i):
                if d[i, j] < 1.15 * (self._rc[atom[i]] + self._rc[atom[j]]):
                    matrix[i, j] = matrix[j, i] = 1

        return matrix

    def get_type(self, atom, matrix):
        n_atom = len(atom)
        htype = [None] * n_atom  # hybridization type
        stype = [0] * n_atom  # sigma profile type
        # 0 : nhb, 1 : OH, 2 : OT, 3 : COOH

        if self.ver == 2002:
            return htype, stype

        for i in range(n_atom):
            ard_i = {j: atom[j] for j in np.flatnonzero(matrix[i])}  # around i
            nabr = ard_i.values()  # neighbor around i
            n_bond = len(nabr)

            if atom[i] == "C":
                stype[i] = 0
                if n_bond == 2:
                    htype[i] = "C(sp)"
                if n_bond == 3:
                    htype[i] = "C(sp2)"
                if n_bond == 4:
                    htype[i] = "C(sp3)"

            if atom[i] == "O":
                stype[i] = 2
                if n_bond == 1:
                    htype[i] = "O(sp2)"
                if n_bond == 2:
                    htype[i] = "O(sp3)"
                    if "H" in nabr:
                        stype[i] = 1

            if atom[i] == "N":
                stype[i] = 2
                if n_bond == 1:
                    htype[i] = "N(sp)"
                if n_bond == 2:
                    htype[i] = "N(sp2)"
                if n_bond >= 3:
                    htype[i] = "N(sp3)"

            if atom[i] == "F":
                stype[i] = 2
                htype[i] = "F"

            if atom[i] == "Cl":
                stype[i] = 0
                htype[i] = "Cl"

            if atom[i] == "H":
                htype[i] = "H(other)"
                if "O" in nabr:
                    stype[i] = 1
                    htype[i] = "H(OH)"

                    # H2O, COOH
                    j = list(ard_i.keys())[0]
                    ard_j = {k: atom[k] for k in np.flatnonzero(matrix[j])}
                    k = [x for x in ard_j.keys() if x != i][0]
                    ard_k = {l: atom[l] for l in np.flatnonzero(matrix[k])}
                    n_O_ard_k = list(ard_k.values()).count("O")

                    if atom[k] == "C" and len(ard_k) == 3 and n_O_ard_k == 2:
                        htype[i] = "H(H2O/COOH)"

                        if self.ver == 2019:
                            stype[i] = 3
                            stype[k] = 3
                            for l in ard_k:
                                if ard_k[l] == "O":
                                    stype[l] = 3

                    if atom[k] == "H":
                        htype[i] = "H(H2O/COOH)"
                        htype[k] = "H(H2O/COOH)"

                if "N" in nabr:
                    stype[i] = 2
                    htype[i] = "H(NH)"

        return htype, stype

    def get_sigma(self, seg, stype):
        n_seg = len(seg)
        segatom, segcoord, segarea, segchrg = (
            seg[:, 0],
            seg[:, 1:4],
            seg[:, 4],
            seg[:, 5],
        )

        d = distance_matrix(segcoord, segcoord)
        r = np.sqrt(segarea / np.pi)

        rcal = r**2 * self._reff**2 / (r**2 + self._reff**2)
        dcal = d**2 / (r**2 + self._reff**2).reshape(-1, 1)

        up = np.einsum("n,n,mn->m", segchrg, rcal, np.exp(-self._fdecay * dcal))
        dn = np.einsum("n,mn->m", rcal, np.exp(-self._fdecay * dcal))
        avechrg = up / dn

        sig = np.linspace(-0.025, 0.025, 51)
        psigA = np.zeros((self._n_psig, 51))

        for i in range(n_seg):
            left = int((avechrg[i] - sig[0]) / 0.001)
            w = (sig[left + 1] - avechrg[i]) / 0.001
            stype_n = stype[int(segatom[i])]
            psigA[stype_n, left] += w * segarea[i]
            psigA[stype_n, left + 1] += (1 - w) * segarea[i]

        if self.ver == 2002:
            return psigA
        else:
            phb = 1 - np.exp(-(sig**2) / 2 / self._sig0**2)
            psigA[0] += np.sum(psigA[1:], axis=0) * (1 - phb)
            psigA[1:] *= phb

            return psigA

    def get_dsp(self, htype):
        if self.ver == 2002 or self.ver == 2010 or None in htype:
            mtype, e = None, None
            return mtype, e

        # molecule type
        n_htype = len(htype)
        mtype = "nhb"
        hba_atoms = ["O(sp2)", "O(sp3)", "N(sp)", "N(sp2)", "N(sp3)", "F"]
        if any(hba_atom in htype for hba_atom in hba_atoms):
            mtype = "hb-a"
        if set(htype) == {"H(other)", "F"}:  # HF
            mtype = "hb-da"
        if "H(OH)" in htype or "H(NH)" in htype:
            mtype = "hb-da"
        if "H(H2O/COOH)" in htype:
            mtype = "COOH"
        if set(htype) == {"H(H2O/COOH)", "O(sp3)", "H(H2O/COOH)"}:  # H2O
            mtype = "H2O"

        # dispersive parameter
        hdict = {
            "C(sp)": 66.0691,
            "C(sp2)": 117.4650,
            "C(sp3)": 115.7023,
            "O(sp2)": -11.0549,
            "O(sp3)": 95.6184,
            "N(sp)": 109.6621,
            "N(sp2)": 84.6268,
            "N(sp3)": 15.4901,
            "F": 52.9318,
            "Cl": 104.2534,
            "H(OH)": 19.3477,
            "H(NH)": 141.1709,
            "H(H2O/COOH)": 58.3301,
            "H(other)": 0,
        }
        e = np.array([hdict[htype[i]] for i in range(n_htype)])
        n = len(e[e != 0])
        e = np.sum(e) / n

        return mtype, e

    def add_comp(self, file, name=None):
        A, V, atom, coord, seg = self.get_cosmo(file)
        matrix = self.get_bond(atom, coord)
        htype, stype = self.get_type(atom, matrix)
        psigA = self.get_sigma(seg, stype)
        mtype, e = self.get_dsp(htype)

        self.A.append(A)
        self.V.append(V)
        self.psigA.append(psigA)
        self.mtype.append(mtype)
        self.e.append(e)
        self.name.append(name)

    def del_comp(self):
        self.A = []
        self.V = []
        self.psigA = []
        self.mtype = []
        self.e = []
        self.name = []

    def ln_gam_comb(self):
        x, A, V = np.array(self.x), np.array(self.A), np.array(self.V)

        q = A / self._q0
        r = V / self._r0
        l = (self._z / 2) * (r - q) - (r - 1)
        tht = q / np.sum(x * q)
        phi = r / np.sum(x * r)

        ln_gam_comb = (
            np.log(phi) + self._z / 2 * q * np.log(tht / phi) + l - phi * np.sum(x * l)
        )

        return ln_gam_comb

    def cal_DelW(self):
        sig = np.linspace(-0.025, 0.025, 51)
        sigT = sig.reshape(-1, 1)

        DelW = np.zeros((self._n_psig, self._n_psig, 51, 51))
        for i in range(self._n_psig):
            for j in range(i + 1):
                if self.ver == 2002:
                    acchb = np.maximum(sig, sigT) - self._sighb
                    donhb = np.minimum(sig, sigT) + self._sighb
                    maxacc = np.where(acchb > 0, acchb, 0)
                    mindon = np.where(donhb < 0, donhb, 0)
                    chb_part = -self._chb[i][j] * maxacc * mindon
                else:
                    mask = (sig * sigT) < 0
                    chb_part = np.where(mask, self._chb[i, j] * (sig - sigT) ** 2, 0)

                DelW[i, j] = DelW[j, i] = (
                    self._cES(self.T) * (sig + sigT) ** 2 - chb_part
                )

        return DelW

    def cal_sigma_mix(self):
        x = np.array(self.x)
        A = np.array(self.A)
        psigA = np.array(self.psigA)

        psigm = np.einsum("i,itm->tm", x, psigA) / np.sum(x * A)
        return psigm

    def ln_gam_res(self):
        A = np.array(self.A)

        psigA = np.array(self.psigA)
        psig = np.einsum("itm,i->itm", psigA, 1 / A)
        psigm = self.cal_sigma_mix()

        DelW = self.cal_DelW()
        exp_DelW = np.exp(-DelW / self._R / self.T)

        Ap = np.einsum("stmn,isn->istmn", exp_DelW, psig)  # A^(+)
        Apm = np.einsum("stmn,sn->stmn", exp_DelW, psigm)  # A^(+)_mix

        for _iter in range(5001):
            if _iter == 5000:
                raise Exception("Convergence Failed")
            if _iter == 0:
                Gam_old = np.ones(np.shape(psig))
                Gamm_old = np.ones(np.shape(psigm))
                Gam = 1 / np.einsum("istmn,isn->itm", Ap, Gam_old)
                Gamm = 1 / np.einsum("stmn,sn->tm", Apm, Gamm_old)
            else:
                Gam_old = Gam
                Gamm_old = Gamm
                Gam = 1 / np.einsum("istmn,isn->itm", Ap, Gam)
                Gamm = 1 / np.einsum("stmn,sn->tm", Apm, Gamm)

                Gam_diff = abs((Gam - Gam_old) / Gam_old)
                Gamm_diff = abs((Gamm - Gamm_old) / Gamm_old)

                Gam = (Gam + Gam_old) / 2
                Gamm = (Gamm + Gamm_old) / 2
                if (np.all(Gam_diff) < 1e-5) and (np.all(Gamm_diff) < 1e-5):
                    break

        Gam_part = np.log(Gamm) - np.log(Gam)
        ln_gam_res = np.einsum("itm,itm->i", psigA, Gam_part) / self._aeff

        return ln_gam_res

    def ln_gam_dsp(self):
        e = np.array(self.e)
        eT = e.reshape(-1, 1)
        n_mole = len(self.x)

        if None in e or None in self.mtype:
            ln_gam_dsp = np.array([None] * n_mole)
            return ln_gam_dsp

        w = np.ones((n_mole, n_mole)) * 0.27027
        wpair = [{"H2O", "hb-a"}, {"COOH", "nhb"}, {"COOH", "hb-da"}, {"H2O", "COOH"}]
        for i in range(n_mole):
            for j in range(i):
                mpair = {self.mtype[i], self.mtype[j]}
                if any(pair == mpair for pair in wpair):
                    w[i][j] = w[j][i] = -0.27027

        A = w * (0.5 * (e + eT) - np.sqrt(e * eT))  # not area

        ln_gam_dsp = np.zeros(n_mole)
        for i in range(n_mole):
            for j in range(n_mole):
                if i != j:
                    ln_gam_dsp[i] += self.x[j] * A[i, j]
                if j > i:
                    ln_gam_dsp[i] -= self.x[i] * self.x[j] * A[i, j]

        return ln_gam_dsp

    def gam(self):
        ln_gam_comb = self.ln_gam_comb()
        ln_gam_res = self.ln_gam_res()
        ln_gam_dsp = self.ln_gam_dsp()

        if None in ln_gam_dsp:
            ln_gam = ln_gam_comb + ln_gam_res
        else:
            ln_gam = ln_gam_comb + ln_gam_res + ln_gam_dsp
        gam = np.exp(ln_gam)

        return gam
