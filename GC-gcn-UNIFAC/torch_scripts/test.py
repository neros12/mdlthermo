from os.path import join as opj
from pathlib import Path

import torch

from get_dataset import retrieve_dataset
from NRTL_layers import VLE_Model


print(opj(Path(__file__).resolve().parents[1], "models", "test1.pt"))


trained_model = torch.load(
    opj(Path(__file__).resolve().parents[1], "models", "test1.pt"),
    map_location=lambda storage, loc: storage,
)

# trained_model = VLE_Model()
# trained_model.load_state_dict(trained_model_dict)

x_data, y_data = retrieve_dataset()
a = x_data[0][8:28]
b = x_data[1][8:28]
c = x_data[2][8:28]
d = x_data[3][8:28]
e = x_data[4][8:28]
f = x_data[5][8:28]
g = x_data[6][8:28]
h = x_data[7][8:28]

result = trained_model(a, b, c, d, e, f, g, h)

pass
