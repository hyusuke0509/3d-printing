import numpy as np, pickle
d=pickle.load(open('mesh.pkl','rb'))['3']; V,T=d['v'],d['t']
st=np.load('tri_state.npy')
def uf_components(n, pairs):
    par=np.arange(n)
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for a,b in pairs:
        ra,rb=find(a),find(b)
        if ra!=rb: par[ra]=rb
    return np.array([find(i) for i in range(n)])
# edge adjacency
E=np.concatenate([T[:,[0,1]],T[:,[1,2]],T[:,[2,0]]]); fid=np.tile(np.arange(len(T)),3)
E.sort(1); key=E[:,0]*10**6+E[:,1]; o=np.argsort(key,kind='stable'); ks=key[o]
same=np.where(ks[1:]==ks[:-1])[0]
pairs=np.stack([fid[o][same],fid[o][same+1]],1)
print('edges shared count, nonmanifold?', np.bincount(np.unique(ks,return_counts=True)[1]))
lab=uf_components(len(T),pairs)
u,c=np.unique(lab,return_counts=True); print('shells',len(u),sorted(c)[::-1][:10])
for r in u[np.argsort(-c)][:10]:
    m=lab==r; vs=V[T[m]].reshape(-1,3)
    print(r,m.sum(),'states',np.bincount(st[m],minlength=3),'bbox',vs.min(0).round(1),vs.max(0).round(1))
w=st!=2
wp=pairs[w[pairs[:,0]]&w[pairs[:,1]]]
lab2=uf_components(len(T),wp)
u2,c2=np.unique(lab2[w],return_counts=True)
print('white regions',len(u2))
for r in u2[np.argsort(-c2)][:25]:
    m=lab2==r; vs=V[T[m]].reshape(-1,3)
    print(r,m.sum(),'shell',np.unique(lab[m]),'bbox',vs.min(0).round(1),vs.max(0).round(1))
np.save('shell.npy',lab); np.save('white.npy',lab2)
