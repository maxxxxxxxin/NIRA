import pandas as pd
import numpy as np
from pyscenic.genesig import GeneSignature
from pyscenic.aucell import create_rankings, enrichment


class Pathstrapping():
    def __init__(self,
                T = 10000,
                binary = False,
                threshold = 0.01):
        self.T = T
        self.threshold = threshold
        self.binary = binary
        
### df is a dataframe which index is the cellnames and coloumn names is genes    
    def fit(self, df, path):
        self.cells = df.index
        self.genes = df.columns
        m, n = df.shape
        self.m = m
        self.n = n

        genes = self.genes      
        rnk_df = create_rankings(df)
        
        gs0 = GeneSignature('default', path)
        auc0 = enrichment(rnk_df, gs0)
        auc0 = auc0['AUC'].values.reshape(-1, 1)
        
        samples = self.bootstrap(path)
        AUC = np.empty((self.m, self.T))
        for idx, sample in enumerate(samples):
            gs = GeneSignature(str(idx), sample)
            auc = enrichment(rnk_df, gs)
            AUC[:,idx] = auc['AUC'].values
            
        self.AUC = AUC
        p = np.sum((AUC - auc0) >= 0, axis=1)/self.T
        P = pd.DataFrame(p)
        P.index = self.cells
        self.Pmat = P
        
        if self.binary == True:
            status = np.where(p < self.threshold, 'On', 'Off')
            status = pd.DataFrame(status)
            status.index = self.cells
            self.binary = status

    
    
    def bootstrap(self, path):
        s = len(path)
        genes = self.genes
        samples = []
        for i in range(self.T):
            sample = np.random.choice(genes, s, replace=False)
            samples.append(sample)
            
        return samples
            
            
            
            
            
            
            
            
            
            
            
        