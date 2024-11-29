import torch
from torch import nn


class GraphConvolution(nn.Module):
    def __init__(self, input_dim=50, output_dim=50):
        super().__init__()
        self.fcn = nn.Sequential(
            nn.Linear(input_dim, output_dim),
            nn.ReLU(),
        )

    def forward(self, nfm, efm):
        x = torch.bmm(efm, nfm)
        x = self.fcn(x)

        return x


class GraphReadOut(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        x = torch.sum(x, 1)

        return x


class VaporPressureModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.gcn1 = GraphConvolution(50, 128)
        self.gcn2 = GraphConvolution(128, 128)
        self.read_out = GraphReadOut()
        self.param_a = nn.Sequential(
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            nn.ReLU(),
        )
        self.param_b = nn.Sequential(
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            nn.ReLU(),
        )
        self.param_c = nn.Sequential(
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            nn.ReLU(),
        )

    def forward(self, nfm, efm, T):
        x = self.gcn1(nfm, efm)
        x = self.gcn2(x, efm)
        x = self.read_out(x)
        a = self.param_a(x)
        b = self.param_b(x)
        c = self.param_c(x)
        lnP = a + b / (T + c)
        P = torch.exp(lnP)

        return P


class NRTLActivityCoefficentModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.gcn1 = GraphConvolution(50, 128)
        self.gcn2 = GraphConvolution(128, 128)
        self.read_out = GraphReadOut()
        self.interaction_parameter = nn.Sequential(
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, 2),
        )

    def NRTL(self, X: torch.Tensor, A: torch.Tensor, B: torch.Tensor, T: torch.Tensor):
        Nc = 2
        alpha = 0.3

        X = X.unsqueeze(1)
        T = T.unsqueeze(2).repeat(1, Nc, Nc)

        tau = A + B / T
        G = torch.exp(-alpha * tau)
        GT = torch.transpose(G, 1, 2)
        tauT = torch.transpose(tau, 1, 2)

        temp0 = torch.bmm(X, tau * G) / torch.bmm(X, G)
        temp1 = torch.transpose(X, 1, 2).repeat(1, 1, Nc) * GT
        temp2 = torch.sum(X.repeat(1, Nc, 1) * GT, 2, keepdim=True).repeat(1, 1, Nc)
        temp3 = torch.sum(X.repeat(1, Nc, 1) * GT * tauT, 2, keepdim=True).repeat(
            1, 1, Nc
        )

        lnGamma = temp0 + torch.sum(
            temp1 / temp2 * (tauT - temp3 / temp2),
            1,
            keepdim=True,
        )
        lnGamma = lnGamma.squeeze(1)

        return lnGamma

    def forward(self, nfm1, efm1, nfm2, efm2, x1, x2, T):
        _batch_size, _ = T.shape
        _Nc = 2

        m1 = self.gcn1(nfm1, efm1)
        m1 = self.gcn2(m1, efm1)
        m1 = self.read_out(m1)

        m2 = self.gcn1(nfm2, efm2)
        m2 = self.gcn2(m2, efm2)
        m2 = self.read_out(m2)

        m12 = torch.cat([m1, m2], 1)
        m12 = self.interaction_parameter(m12)
        m21 = torch.cat([m2, m1], 1)
        m21 = self.interaction_parameter(m21)

        a_matrix = torch.zeros(_batch_size, _Nc, _Nc)
        a_matrix[:, 0, 1] = m12[:, 0]
        a_matrix[:, 1, 0] = m21[:, 0]

        b_matrix = torch.zeros(_batch_size, _Nc, _Nc)
        b_matrix[:, 0, 1] = m12[:, 1]
        b_matrix[:, 1, 0] = m21[:, 1]

        x = torch.cat([x1, x2], 1)
        lnGamma = self.NRTL(x, a_matrix, b_matrix, T)
        gamma = torch.exp(lnGamma)

        return gamma


class VLE_Model(nn.Module):
    def __init__(self):
        super().__init__()
        self.activity = NRTLActivityCoefficentModel()
        self.Psat = VaporPressureModel()

    def forward(self, nfm1, efm1, nfm2, efm2, x1, x2, T, P):
        activity_result = self.activity(nfm1, efm1, nfm2, efm2, x1, x2, T)
        Psat1 = self.Psat(nfm1, efm1, T)
        Psat2 = self.Psat(nfm2, efm2, T)

        # P_estimated = (
        #     x1 * activity_result[:, 0].unsqueeze(-1) * Psat1
        #     + x2 * activity_result[:, 1].unsqueeze(-1) * Psat2
        # )
        y1_estimated = x1 * activity_result[:, 0].unsqueeze(-1) * Psat1 / P
        y2_estimated = x2 * activity_result[:, 1].unsqueeze(-1) * Psat2 / P

        return torch.cat([y1_estimated, y2_estimated], 1)


class CustomLoss(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, y_pred, y_exp):
        y_temp_exp = y_exp[:, 0:2]
        loss = torch.sum(((y_pred - y_temp_exp) / y_temp_exp) ** 2)

        return loss


# PyTorch Init
if __name__ == "__main__":
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        device_name = torch.cuda.get_device_name(torch.cuda.current_device())
    else:
        device_name = "cpu"
    torch.set_default_device(device)
