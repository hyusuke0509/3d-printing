import numpy as np, pickle, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
d=pickle.load(open('mesh.pkl','rb'))['3']; V,T=d['v'],d['t']
wl=np.load('white.npy'); m=wl==316163
C=np.load('eggfit.npy'); R=C[3]; C=C[:3]
c=V[T].mean(1); dall=np.linalg.norm(c-C,axis=1)-R
E=np.concatenate([T[:,[0,1]],T[:,[1,2]],T[:,[2,0]]]); fid=np.tile(np.arange(len(T)),3)
E.sort(1); key=E[:,0]*10**6+E[:,1]; o=np.argsort(key,kind='stable'); ks=key[o]
same=np.where(ks[1:]==ks[:-1])[0]; pairs=np.stack([fid[o][same],fid[o][same+1]],1)
np.save('pairs.npy',pairs)
import sys
th=float(sys.argv[1]) if len(sys.argv)>1 else 0.3
cand=m&(dall>th)&(c[:,2]>-26)
pp=pairs[cand[pairs[:,0]]&cand[pairs[:,1]]]
G=sp.coo_matrix((np.ones(len(pp)),(pp[:,0],pp[:,1])),shape=(len(T),len(T)))
k,lab=connected_components(G,directed=False)
u,cnt=np.unique(lab[cand],return_counts=True)
for r in u[np.argsort(-cnt)][:6]:
    s=lab==r; print(r,s.sum(),c[s].min(0).round(1),c[s].max(0).round(1))
best=u[np.argmax(cnt)]
np.save('lwing_mask.npy',lab==best)
