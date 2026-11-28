import numpy as np, pickle, collections
def decode(s, V):
    """V: 3x3 array of triangle vertex coords. returns list of (state, 3x3 coords)"""
    nib=[int(c,16) for c in reversed(s)]
    pos=[0]
    out=[]
    def nxt():
        n=nib[pos[0]]; pos[0]+=1; return n
    def rec(v):
        code=nxt(); ns=code&3; 
        if ns==0:
            st=code>>2
            if st==3: st=nxt()+3
            out.append((st,v)); return
        sp=code>>2
        a=[v[(sp+j)%3] for j in range(3)]
        if ns==1:
            m=(a[1]+a[2])/2; vv=[a[0],a[1],m,a[2]]
            ch=[(vv[0],vv[1],vv[2]),(vv[2],vv[3],vv[0])]
        elif ns==2:
            m01=(a[0]+a[1])/2; m02=(a[0]+a[2])/2
            vv=[a[0],m01,a[1],a[2],m02]
            ch=[(vv[0],vv[1],vv[4]),(vv[1],vv[2],vv[4]),(vv[2],vv[3],vv[4])]
        else:
            m01=(a[0]+a[1])/2; m12=(a[1]+a[2])/2; m20=(a[2]+a[0])/2
            vv=[a[0],m01,a[1],m12,a[2],m20]
            ch=[(vv[0],vv[1],vv[5]),(vv[1],vv[2],vv[3]),(vv[3],vv[4],vv[5]),(vv[1],vv[3],vv[5])]
        res=[None]*len(ch)
        for i in reversed(range(len(ch))):
            rec(np.array(ch[i]))
    rec(np.array(V))
    return out
if __name__=='__main__':
    d=pickle.load(open('mesh.pkl','rb'))['3']
    V,T,P=d['v'],d['t'],d['p']
    cnt=collections.Counter(); area=collections.Counter()
    tri_state=np.zeros(len(T),int)  # dominant state
    mixed=0
    for i,p in enumerate(P):
        if not p: continue
        leaves=decode(p,V[T[i]])
        a=collections.Counter()
        for st,v in leaves:
            ar=np.linalg.norm(np.cross(v[1]-v[0],v[2]-v[0]))/2
            a[st]+=ar; area[st]+=ar
        if len(a)>1: mixed+=1
        tri_state[i]=a.most_common(1)[0][0]
    print('area by state',area,'mixed tris',mixed)
    print(collections.Counter(tri_state.tolist()))
    np.save('tri_state.npy',tri_state)
