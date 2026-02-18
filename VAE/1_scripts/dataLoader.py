# define a Dataset object for data loading and train / test splitting
import random
import math
import torch
from torch.utils.data import Dataset, DataLoader

class GWASDataset(Dataset):
    def __init__(self, labelfile, datafile, train=True, train_split=0.8, random_seed=1):
        self.train_split = train_split
        self.random_seed = random_seed
        
        self.data = load_dataset(labelfile, datafile, train, train_split, random_seed)
    
    def __len__(self):
        return self.data.shape[0]
    
    def __getitem__(self, idx):
        #l = self.labels[idx]
        d = self.data[idx]
        return d#, l
    
def load_dataset(labelfile, datafile, train, train_split, random_seed):
    all_data = torch.load(datafile)
    n = all_data.shape[0]
    
    random.seed(random_seed)
    train_idxs = random.sample(range(n), math.floor(n*train_split))
    test_idxs = [i for i in range(n) if i not in train_idxs]
    if train:
        idxs=train_idxs
    else:
        idxs=test_idxs
    
    return all_data[idxs, :]

def setup_data_loaders(batch_size=128, use_cuda=False, data_dir = "../0_data/", datafile="data.pt", labelfile="labels.csv", random_seed=314):
    train_set = GWASDataset(data_dir+labelfile, data_dir+datafile,
                            train=True, random_seed=random_seed)
    test_set = GWASDataset(data_dir+labelfile, data_dir+datafile,
                           train=False, random_seed=random_seed)
    
    kwargs = {"num_workers":0, "pin_memory": use_cuda}
    train_loader = DataLoader(dataset=train_set,
                              batch_size=batch_size, shuffle=True, **kwargs)
    test_loader = DataLoader(dataset=test_set,
                             batch_size=batch_size, shuffle=False, **kwargs)
    
    return train_loader, test_loader
    
        
        
        
        
        