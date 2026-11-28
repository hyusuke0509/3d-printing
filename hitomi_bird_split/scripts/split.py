import numpy as np, pickle, collections, scipy.sparse as sp, triangle, shapely.geometry as sg
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
import manifold3d as mf
d=pickle.load(open('mesh.pkl','rb'))['3']; V,T=d['v'],d['t']
wl=np.load('white.npy'); pairs=np.load('pairs.npy'); nT=len(T)
A=sp.coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(nT,nT)); A=(A+A.T).tocsr()
masks={'tail':wl==340019,'wing_R':wl==322110,'wing_L':np.load('lwing_mask.npy')}
def comps(mask):
    idx=np.where(mask)[0]; sub=A[idx][:,idx]; k,l=connected_components(sub,directed=False)
    return idx,l,k
def clean(mask):
    mask=mask.copy()
    # morphological closing then opening (1 ring)
    for _ in range(2): mask|= (A@mask.astype(float))>=2
    for _ in range(2): mask&= (A@mask.astype(float))>=2
    # keep largest comp
    idx,l,k=comps(mask); big=np.bincount(l).argmax(); mask[:]=False; mask[idx[l==big]]=True
    # fill holes: complement comps other than the largest
    idx,l,k=comps(~mask); big=np.bincount(l).argmax()
    # complement may have small separate shell (968) -> only fill comps adjacent to mask
    for j in range(k):
        if j==big: continue
        t=idx[l==j]
        if (A[t]@mask.astype(float)).sum()>0: mask[t]=True
    return mask
def boundary_loops(F):
    E=np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]])
    s=set(map(tuple,E.tolist())); b=[e for e in map(tuple,E.tolist()) if (e[1],e[0]) not in s]
    nxt=collections.defaultdict(list)
    for a,c in b: nxt[a].append(c)
    pinch=[v for v,l in nxt.items() if len(l)>1]
    used=set(); L=[]
    for a,c in b:
        if (a,c) in used: continue
        loop=[a]; cur=(a,c)
        while cur not in used:
            used.add(cur); loop.append(cur[1]); cands=[x for x in nxt[cur[1]] if (cur[1],x) not in used]
            if not cands: break
            cur=(cur[1],cands[0])
        L.append(loop[:-1])
    return L,pinch
def make_cap(loop):
    P=V[loop]; o=P.mean(0); u,s,vt=np.linalg.svd(P-o)
    rng=np.random.default_rng(0); dirs=np.r_[vt[2:3],vt[2:3]]; dirs/=np.linalg.norm(dirs,axis=1,keepdims=True)
    dirs=dirs[np.argsort(-np.abs(dirs@vt[2]))]
    for nrm in dirs:
        e1=np.cross(nrm,[0.3,0.5,0.8]); e1/=np.linalg.norm(e1); e2=np.cross(nrm,e1)
        Q=np.c_[(P-o)@e1,(P-o)@e2]
        if sg.Polygon(Q).is_valid: break
        nrm=None; break
    if nrm is None:
        # map loop to circle by arc length (disk domain)
        seglen=np.linalg.norm(P-np.roll(P,-1,0),axis=1); t=np.r_[0,np.cumsum(seglen)[:-1]]/seglen.sum()*2*np.pi
        rad=seglen.sum()/(2*np.pi); Q=np.c_[rad*np.cos(t),rad*np.sin(t)]; e1=vt[0]; e2=vt[1]
        print('   disk mapping')
    else: print('   proj normal',nrm.round(2),'dot pca',abs(nrm@vt[2]).round(2))
    n=len(loop); seg=np.c_[np.arange(n),(np.arange(n)+1)%n]
    edge=np.mean(np.linalg.norm(np.diff(np.r_[Q,Q[:1]],axis=0),axis=1))
    tr=triangle.triangulate({'vertices':Q,'segments':seg},f'pq20Ya{edge**2*3:.6f}')
    TV=tr['vertices']; TT=tr['triangles']; assert len(TV)>=n and np.allclose(TV[:n],Q)
    m=len(TV); X=np.zeros((m,3)); X[:n]=P; X[n:]=o
    # harmonic membrane for interior verts
    if m>n:
        E=np.r_[TT[:,[0,1]],TT[:,[1,2]],TT[:,[2,0]]]
        L=sp.coo_matrix((np.ones(len(E)),(E[:,0],E[:,1])),shape=(m,m)); L=((L+L.T)>0).astype(float).tocsr()
        deg=np.asarray(L.sum(1)).ravel(); Lap=sp.diags(deg)-L
        I=np.arange(n,m); B=np.arange(n)
        for k in range(3):
            X[I,k]=spsolve(Lap[I][:,I].tocsc(),-(Lap[I][:,B]@X[B,k]))
    gidx=np.r_[np.array(loop),-1-np.arange(m-n)]  # negative = new verts
    return X[n:],gidx[TT]
def assemble(tri_list_old, cap_tris, newV, flip):
    # build mesh with vertices: old used + new
    ct=cap_tris[:,::-1] if flip else cap_tris
    allt=np.r_[tri_list_old,ct]
    used=np.unique(allt[allt>=0]); remap={}
    oldmap=-np.ones(len(V),int); oldmap[used]=np.arange(len(used))
    VV=np.r_[V[used],newV]
    tt=np.where(allt>=0,oldmap[np.clip(allt,0,None)],len(used)+(-1-allt))
    return VV,tt
def vol(VV,tt):
    a,b,c=VV[tt[:,0]],VV[tt[:,1]],VV[tt[:,2]]; return (a*np.cross(b,c)).sum()/6
def M(VV,tt): return mf.Manifold(mf.Mesh(vert_properties=VV.astype(np.float32),tri_verts=tt.astype(np.uint32)))
orig=M(V,T); print('orig status',orig.status(),'vol',orig.volume(),vol(V,T))
rest=np.ones(nT,bool); caps=[]; parts={}
for name,mask in masks.items():
    mask=clean(mask)&rest; rest&=~mask
    L,pinch=boundary_loops(T[mask]); print(name,'tris',mask.sum(),'loops',[len(l) for l in L],'pinch',len(pinch))
    assert len(L)==1 and not pinch
    newV,ct=make_cap(L[0]); print('   cap tris',len(ct))
    # feather: region + cap ; cap tris index loop verts in same order as region boundary => region boundary edges a->b ; cap must have b->a
    l=L[0]; e0=(l[0],l[1])
    has=any(((t[0],t[1])==e0 or (t[1],t[2])==e0 or (t[2],t[0])==e0) for t in ct.tolist())
    flip=has
    VV,tt=assemble(T[mask],ct,newV,flip)
    parts[name]=(VV,tt); caps.append((ct,newV,flip))
    f=M(VV,tt); print('  feather status',f.status(),'vol',round(vol(VV,tt),2),'∩orig',round((f^orig).volume(),2))
# body: rest + all caps reversed
allT=[T[rest]]; allV=[V]; off=len(V)
for ct,newV,flip in caps:
    ct2=ct[:,::-1] if not flip else ct
    ct2=np.where(ct2>=0,ct2,off+(-1-ct2)); allT.append(ct2); allV.append(newV); off+=len(newV)
BV=np.concatenate(allV); BT=np.concatenate(allT)
used=np.unique(BT); mp=-np.ones(len(BV),int); mp[used]=np.arange(len(used)); BV=BV[used]; BT=mp[BT]
b=M(BV,BT); print('body status',b.status(),'vol',round(vol(BV,BT),2),'∩orig',round((b^orig).volume(),2))
for n,(VV,tt) in parts.items():
    print(n,'∩body',round((M(VV,tt)^b).volume(),3))
parts['body']=(BV,BT); np.save('rest.npy',rest)
pickle.dump(parts,open('parts.pkl','wb'))
