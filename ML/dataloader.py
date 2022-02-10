# %%
from torch.utils.data import Dataset
from torchvision.io import read_image
import matplotlib.pyplot as plt
import os 
import pandas as pd 

# %% 
class StarImageDataset(Dataset):
    def __init__(self, data_file, img_dir, transform=None, target_transform=None):
        self.img_labels = pd.read_csv(data_file)
        self.img_dir = img_dir
        self.transform = transform
        self.target_transform = target_transform

    def __len__(self):
        return len(self.img_labels)
    
    def __getitem__(self, idx):
        img_path = os.path.join(self.img_dir, self.img_labels.iloc[idx, 0])
        image = read_image(img_path)
        label = self.img_labels.iloc[idx, 1]
        if self.transform: 
            image = self.transform(image)
        if self.target_transform:
            label = self.target_transform(label)
        return image, label 
