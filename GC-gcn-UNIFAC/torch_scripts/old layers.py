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
        x = torch.matmul(efm, nfm)
        x = self.fcn(x)

        return x


class GraphReadOut(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        x = torch.sum(x, 1)

        return x


class RLayer(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.diagonal_weights = nn.Parameter(torch.randn(dim) + 10)

    def forward(self, x):
        diagonal_matrix = torch.diag(self.diagonal_weights)
        x = torch.matmul(x, diagonal_matrix)

        return x


class QLayer(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.diagonal_weights = nn.Parameter(torch.randn(dim) + 10)

    def forward(self, x):
        diagonal_matrix = torch.diag(self.diagonal_weights)
        x = torch.matmul(x, diagonal_matrix)

        return x


class InteractionLayer(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.weights = nn.Parameter((torch.randn(dim, dim) * 0.01).fill_diagonal_(0))
        self.dim = dim

    def forward(self):
        interaction_parameters = self.weights * (1 - torch.eye(self.dim))

        return interaction_parameters


class InteractionLayer(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.fcn = nn.Sequential(
            nn.Linear(dim, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
        )

    def forward(self, x):
        x = self.fcn(x)

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
        )
        self.param_b = nn.Sequential(
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
        )
        self.param_c = nn.Sequential(
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
        )

    def forward(self, nfm, efm, T):
        x = self.gcn1(nfm, efm)
        x = self.gcn2(x, efm)
        x = self.read_out(x)
        a = self.param_a(x)
        b = self.param_b(x)
        c = self.param_c(x)
        lnP = a - b / (T + c)
        P = torch.exp(lnP)

        return P


class ActivityCoefficentModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.gcn1 = GraphConvolution(50, 128)
        self.gcn2 = GraphConvolution(128, 128)
        self.read_out = GraphReadOut()
        self.r_layer = RLayer(128)
        self.q_layer = QLayer(128)
        self.inter1 = InteractionLayer(128)

    def forward(self, nfm1, efm1, nfm2, efm2, x1, x2, T):
        _batch_size, _ = T.shape
        dim = 128

        p1 = self.gcn1(nfm1, efm1)
        p1 = self.gcn2(p1, efm1)
        p1 = self.read_out(p1)
        r1 = self.r_layer(p1)
        r1 = torch.sum(r1, 1).unsqueeze(-1)
        q1 = self.q_layer(p1)
        q1 = torch.sum(q1, 1).unsqueeze(-1)
        e1 = self.r_layer(p1) / q1

        p2 = self.gcn1(nfm2, efm2)
        p2 = self.gcn2(p2, efm2)
        p2 = self.read_out(p2)
        r2 = self.r_layer(p2)
        r2 = torch.sum(r2, 1).unsqueeze(-1)
        q2 = self.q_layer(p2)
        q2 = torch.sum(q2, 1).unsqueeze(-1)
        e2 = self.r_layer(p2) / q2

        J1 = r1 / (r1 * x1 + r2 * x2)
        J2 = r2 / (r1 * x1 + r2 * x2)
        L1 = q1 / (q1 * x1 + q2 * x2)
        L2 = q2 / (q1 * x1 + q2 * x2)
        theta = (x1 * q1 * e1 + x2 * q2 * e2) / (x1 * q1 + x2 * q2)

        beta1 = torch.exp(self.inter1(e1) / T)
        beta2 = torch.exp(self.inter1(e2) / T)
        s = torch.exp(self.inter1(theta) / T)
        # tow = torch.exp(
        #     -self.inter1().unsqueeze(0).expand(_batch_size, dim, dim)
        #     / T.unsqueeze(-1).expand(_batch_size, dim, dim)
        # )
        # beta1 = torch.matmul(e1.view(-1, 1, dim), tow).view(-1, dim)
        # beta2 = torch.matmul(e2.view(-1, 1, dim), tow).view(-1, dim)
        # s = torch.matmul(theta.view(-1, 1, dim), tow).view(-1, dim)

        ln_r1_R = q1 * (
            1
            - torch.sum((theta * beta1 / s - e1 * torch.log(beta1 / s)), 1).unsqueeze(
                -1
            )
        )
        ln_r2_R = q2 * (
            1
            - torch.sum((theta * beta2 / s - e2 * torch.log(beta2 / s)), 1).unsqueeze(
                -1
            )
        )
        ln_r1_C = 1 - J1 + torch.log(J1) - 5 * q1 * (1 - J1 / L1 + torch.log(J1 / L1))
        ln_r2_C = 1 - J2 + torch.log(J2) - 5 * q2 * (1 - J2 / L2 + torch.log(J2 / L2))

        ln_r1 = ln_r1_R + ln_r1_C
        ln_r2 = ln_r2_R + ln_r2_C

        return torch.cat([torch.exp(ln_r1), torch.exp(ln_r2)], 1)


class VLE_Model(nn.Module):
    def __init__(self):
        super().__init__()
        self.activity = ActivityCoefficentModel()
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
        y_temp_exp = y_exp[:, 0:1]
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
