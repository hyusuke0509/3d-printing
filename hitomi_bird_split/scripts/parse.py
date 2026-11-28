import sys, re, numpy as np, pickle
from lxml import etree
src=sys.argv[1]
ns='{http://schemas.microsoft.com/3dmanufacturing/core/2015/02}'
objs={}
cur=None
for ev,el in etree.iterparse(src,events=('start','end')):
    tag=el.tag.replace(ns,'')
    if ev=='start' and tag=='object':
        cur=el.get('id'); objs[cur]={'v':[],'t':[],'p':[]}
    elif ev=='end' and tag=='vertex':
        objs[cur]['v'].append((float(el.get('x')),float(el.get('y')),float(el.get('z')))); el.clear()
    elif ev=='end' and tag=='triangle':
        objs[cur]['t'].append((int(el.get('v1')),int(el.get('v2')),int(el.get('v3'))))
        objs[cur]['p'].append(el.get('paint_color') or ''); el.clear()
out={k:{'v':np.array(o['v']),'t':np.array(o['t'],dtype=np.int64),'p':o['p']} for k,o in objs.items() if o['v']}
pickle.dump(out,open('mesh.pkl','wb'))
for k,o in out.items(): print(k,o['v'].shape,o['t'].shape,sum(1 for p in o['p'] if p))
