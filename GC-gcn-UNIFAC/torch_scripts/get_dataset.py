import torch


from utils import retrieve_data


def retrieve_dataset():
    x_data, y_data = retrieve_data()
    nfm1, efm1, nfm2, efm2, list_T, list_P, list_x1 = x_data
    nfm1_tensor = torch.tensor(nfm1, dtype=torch.float32)
    efm1_tensor = torch.tensor(efm1, dtype=torch.float32)
    nfm2_tensor = torch.tensor(nfm2, dtype=torch.float32)
    efm2_tensor = torch.tensor(efm2, dtype=torch.float32)
    T_tensor = torch.tensor(list_T, dtype=torch.float32).view(-1, 1)
    P_tensor = torch.tensor(list_P, dtype=torch.float32).view(-1, 1)
    x1_tensor = torch.tensor(list_x1, dtype=torch.float32).view(-1, 1)
    x2_tensor = 1 - x1_tensor
    y1_tensor = torch.tensor(y_data, dtype=torch.float32).view(-1, 1)
    y2_tensor = 1 - y1_tensor

    x_tensor_data = [
        nfm1_tensor,
        efm1_tensor,
        nfm2_tensor,
        efm2_tensor,
        x1_tensor,
        x2_tensor,
        T_tensor,
        P_tensor,
    ]

    y_tensor_data = [y1_tensor, y2_tensor, P_tensor]

    return x_tensor_data, y_tensor_data


# PyTorch Init
if __name__ == "__main__":
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        device_name = torch.cuda.get_device_name(torch.cuda.current_device())
    else:
        device_name = "cpu"
    torch.set_default_device(device)
