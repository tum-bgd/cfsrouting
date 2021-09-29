import numpy as np
import cfsrouting as cfs

class cfg:
    size=25
    n=5


if __name__=="__main__":
    m = np.zeros(cfg.size*cfg.size).reshape(-1,cfg.size).astype(np.uint8);
    p = np.random.uniform(low=0, high=cfg.size, size=cfg.n*2).reshape(-1,2).astype(np.int32)
    m[cfg.size//2,cfg.size//2] = 1
    print(cfs.empty(m,p,0))
    handle = cfs.graphFromBitmap(m,0)
    cfs.dijkstra(handle, [10,10])
    print(cfs.getpath(handle, [15,15],2))
    cfs.dijkstra(handle, [10,10])
    print(cfs.getpath(handle, [15,15],2))
    
    
    
    
