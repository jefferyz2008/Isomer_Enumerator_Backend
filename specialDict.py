from collections import Counter
class dictWithCounts:
    def __init__(self,initial=None):
        self.dict=dict(initial or {})
        self.counts=Counter(self.dict.values())

    def __setitem__(self,k,v):
        old=self.dict.get(k)
        if old is not None:
            self.counts[old]-=1
            if self.counts[old]==0:
                del self.counts[old]
        self.dict[k]=v
        self.counts[v]+=1
    def hasValue(self,v):
        return self.counts.get(v,0)
    def __getitem__(self, k):      
        return self.dict[k]