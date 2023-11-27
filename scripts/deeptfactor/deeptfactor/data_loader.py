import numpy as np
from torch.utils.data import Dataset


class EnzymeDataset(Dataset):
    def __init__(self, data_X, data_Y):
        self.getAAmap()
        self.data_X = data_X
        self.data_Y = data_Y
        
    
    def __len__(self):
        return len(self.data_X)
    
    
    def getAAmap(self):
        aa_vocab = ['A', 'C', 'D', 'E', 
                    'F', 'G', 'H', 'I', 
                    'K', 'L', 'M', 'N', 
                    'P', 'Q', 'R', 'S',
                    'T', 'V', 'W', 'X', 
                    'Y', '_']
        map = {}
        for i, char in enumerate(aa_vocab):
            baseArray = np.zeros(len(aa_vocab)-1)
            if char != '_':
                baseArray[i] = 1
            map[char] = baseArray
        self.map = map
        return

        
    def convert2onehot(self, single_seq):
        single_onehot = []
        for x in single_seq:
            single_onehot.append(self.map[x])
        return np.asarray(single_onehot)
    
    
    def __getitem__(self, idx):
        x = self.data_X[idx]
        x = self.convert2onehot(x)
        y = self.data_Y[idx]
        
        x = x.reshape((1,) + x.shape)
        y = y.reshape(-1)
        return x, y