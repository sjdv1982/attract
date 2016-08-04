import numpy as np
import sys, os



if '--merge' in sys.argv:
  pregrid1=sys.argv[1]
  pregrid2=sys.argv[2]
  
  grid1=np.load(pregrid1)
  grid2=np.load(pregrid2)
  
  if np.shape(grid1)[-2]==np.shape(grid2)[-2]:
    if len(np.shape(grid1))<len(np.shape(grid2)):
      grid1=np.reshape(grid1,grid2.shape[:-1]+(grid1.shape[-1],))
    elif len(np.shape(grid2))<len(np.shape(grid1)):
      grid2=np.reshape(grid2,grid1.shape[:-1]+(grid2.shape[-1],))

    newgrid=np.append(grid1,grid2,axis=-1)
    if '--output' in sys.argv:
      name=sys.argv[sys.argv.index('--output')+1]
    else:
      name='pregrid_merge_'+os.path.splitext(pregrid1)[0]+os.path.splitext(pregrid2)[0]
  else:
    print 'grids with shapes ', np.shape(grid1), ' and ', np.shape(grid2), ' do not fit together'
    sys.exit()
  
else:
  pregrid=sys.argv[1]
  steps=sys.argv[2]
  name=sys.argv[3]

  r=[0.]
  for i in range(4,int(steps)+4):
    r.append(float(sys.argv[i]))
    name+='-'+sys.argv[i]
    
  Pregrid=np.load(pregrid)
  strux=len(Pregrid[0])
  attyp=np.shape(Pregrid)[-1]
  newgrid=np.zeros((int(steps),strux,attyp),dtype=np.uint32)

  for i in range(int(steps)):
    rmin=r[i]
    rmax=r[i+1]
    newgrid[i]=np.sum(Pregrid[rmin:rmax,:,:],axis=0)
    
np.save(name+'.npy', newgrid)
