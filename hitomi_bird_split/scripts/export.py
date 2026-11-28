import numpy as np, pickle, os, struct
out=os.path.expanduser('/home/user/3d-printing/hitomi_bird_split'); os.makedirs(out,exist_ok=True)
parts=pickle.load(open('parts.pkl','rb')); rest=np.load('rest.npy'); st=np.load('tri_state.npy')
dz=np.array([0,0,36.0])
names={'body':'body_本体','tail':'feather_tail_尾羽','wing_R':'feather_wing_R_右翼','wing_L':'feather_wing_L_左翼'}
def stl(path,VV,tt):
    tri=VV[tt]; n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]); n/=np.linalg.norm(n,axis=1,keepdims=True)+1e-12
    rec=np.zeros(len(tt),dtype=[('n','<f4',3),('v','<f4',(3,3)),('a','<u2')]); rec['n']=n; rec['v']=tri
    with open(path,'wb') as f: f.write(b'hitomi bird split'.ljust(80,b' ')); f.write(struct.pack('<I',len(tt))); f.write(rec.tobytes())
with open(f'{out}/hitomi_bird_split.mtl','w') as f:
    for m,c in [('white',(0.95,0.95,0.95)),('black',(0.05,0.05,0.05)),('cut_face',(0.9,0.2,0.2))]:
        f.write(f'newmtl {m}\nKd {c[0]} {c[1]} {c[2]}\n\n')
lines=['mtllib hitomi_bird_split.mtl']; off=1
for k in ['body','tail','wing_R','wing_L']:
    VV,tt=parts[k]; VV=VV+dz
    stl(f'{out}/{k}.stl',VV,tt)
    nb=rest.sum()
    if k=='body': mat=np.r_[np.where(st[rest]==2,'black','white'),np.full(len(tt)-nb,'cut_face')]
    else:
        nr=len(tt)-0; mat=np.full(len(tt),'white')
    lines.append(f'o {k}')
    lines += [f'v {x:.5f} {y:.5f} {z:.5f}' for x,y,z in VV]
    cur=None
    for t,m in zip(tt+off,mat):
        if m!=cur: lines.append(f'usemtl {m}'); cur=m
        lines.append(f'f {t[0]} {t[1]} {t[2]}')
    off+=len(VV)
open(f'{out}/hitomi_bird_split.obj','w').write('\n'.join(lines)+'\n')
os.system(f'ls -la {out}')
