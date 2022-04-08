# %%
import os 
import torch 
from torch import nn 
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
import torch.optim as optim
from dataloader import StarImageDataset
from constants import *
from tqdm import tqdm
from torch.utils.tensorboard import SummaryWriter
# %%
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'Using {device} device')

# %%
class StarNeuralNetwork(nn.Module):
    def __init__(self):
        super(StarNeuralNetwork, self).__init__()
        self.stack = nn.Sequential(
            nn.Conv2d(3, 6, 5),
            nn.ReLU(),
            nn.MaxPool2d(2, 2),
            nn.Conv2d(6, 16, 5),
            nn.ReLU(),
            nn.MaxPool2d(2, 2),
            nn.Conv2d(16, 16, 5),
            nn.ReLU(),
            nn.MaxPool2d(2, 2),
            nn.Flatten(),
            nn.Linear(16 * 121 * 121, 120),
            nn.ReLU(),
            nn.Linear(120, 84),
            nn.ReLU(),
            nn.Linear(84, 10),
            nn.ReLU(),
            nn.Linear(10, 3)
        )

    def forward(self, x):
        return self.stack(x)
net = StarNeuralNetwork()
if device == 'cuda':
    net.to(device)
# %%
star_dataset = StarImageDataset(PATH_TO_DATA_CSV, DATA_PATH)
dataloader = DataLoader(star_dataset, batch_size=BATCH_SIZE, shuffle=True)
# %%
criterion = nn.MSELoss()
optimizer = optim.SGD(net.parameters(), lr=LEARNING_RATE, momentum=MOMENTUM)

#%%
writer = SummaryWriter()
for epoch in tqdm(range(NUM_EPOCHS), position = 0, desc = "Epochs"):
    for i, data in enumerate(tqdm(dataloader, position = 1, desc=f'Batches for Epoch {epoch}', leave = False)):
        inputs, labels = data
        if device == 'cuda':
            inputs, labels = inputs.to(device), labels.to(device)
        optimizer.zero_grad()
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        writer.add_scalar('Loss/train', loss.item(), epoch * len(dataloader) + i)