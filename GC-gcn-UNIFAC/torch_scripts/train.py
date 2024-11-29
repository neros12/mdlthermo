from typing import Optional

import numpy as np
import torch

from get_dataset import retrieve_dataset
from NRTL_layers import CustomLoss, VLE_Model


# PyTorch Init
if __name__ == "__main__":
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        device_name = torch.cuda.get_device_name(torch.cuda.current_device())
    else:
        device_name = "cpu"
    torch.set_default_device(device)


def mini_batching(x_data, y_data, batch_size):
    if len(x_data[0]) != len(y_data[0]):
        raise ValueError("Size of x_data and y_data are not equal!")

    data_length = len(x_data[0])
    mini_batches = []
    for i in range(0, data_length, batch_size):
        x_batch = [
            x_data[0][i : i + batch_size],
            x_data[1][i : i + batch_size],
            x_data[2][i : i + batch_size],
            x_data[3][i : i + batch_size],
            x_data[4][i : i + batch_size],
            x_data[5][i : i + batch_size],
            x_data[6][i : i + batch_size],
            x_data[7][i : i + batch_size],
        ]
        y_batch = torch.cat(
            [
                y_data[0][i : i + batch_size],
                y_data[1][i : i + batch_size],
                y_data[2][i : i + batch_size],
            ],
            1,
        )
        mini_batches.append([x_batch, y_batch])

    return mini_batches


def train_loop(data_loader, ML_model, optimizer, verbose: bool):
    losses = []
    loss_fn = CustomLoss()
    ML_model.train()
    for batch, (x_data, y_data) in enumerate(data_loader):

        pred_y = ML_model(*x_data)
        loss = loss_fn(pred_y, y_data[:, 0:2])
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()

        losses.append(loss.item())
        # if verbose:
        #     print(f"Batch {batch+1}/{len(data_loader)}, Loss: {loss.item():.4f}")

    return np.mean(losses)


def train_model(
    batch_size: int,
    epochs: int,
    seed: Optional[int] = None,
    verbose=True,
):
    x_data, y_data = retrieve_dataset()
    data_loader = mini_batching(x_data, y_data, batch_size)
    ML_model = VLE_Model()

    if verbose:
        print("-------------------------------------------")
        print("             PYTORCH ACTIVATED             ")
        print(f"at: {device} ({device_name})")
        print("-------------------------------------------")
        print("")

    optimizer = torch.optim.Adam(ML_model.parameters())
    train_losses = []
    for epoch in range(epochs):
        if verbose:
            print(f"Epoch {epoch+1}/{epochs}")

        train_loss = train_loop(data_loader, ML_model, optimizer, verbose)
        print("loss :", train_loss)
        train_losses.append(train_loss)

    return ML_model


model = train_model(batch_size=16384, epochs=3000)
torch.save(model, "test.pt")
pass
