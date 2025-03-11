import numpy as np
import h5py

if __name__ == '__main__':
    
    D_id = 'AlbertYoung_chorus'
    
    with h5py.File(D_id + '.h5', 'w') as f:
        
        nalpha0 = 90
        nE = 49
        nL = 40
        
        alpha0 = np.linspace(1, 90, nalpha0)
        logE = np.linspace(np.log(0.1), np.log(5.0), nE)
        L = np.linspace(3, 8, nL)
        
        f.create_dataset('alpha0', data=alpha0)
        f.create_dataset('E', data=np.exp(logE))
        f.create_dataset('L', data=L)
        
        daa2d = np.loadtxt(D_id + "/" + D_id + ".Daa")
        dap2d = np.loadtxt(D_id + "/" + D_id + ".Dap")
        dpp2d = np.loadtxt(D_id + "/" + D_id + ".Dpp")
        
        daa = np.repeat(daa2d[:, :, np.newaxis], nL, axis=2)
        dap = np.repeat(dap2d[:, :, np.newaxis], nL, axis=2)
        dpp = np.repeat(dpp2d[:, :, np.newaxis], nL, axis=2)
        
        f.create_dataset('Daa', data=daa)
        f.create_dataset('Dap', data=dap)
        f.create_dataset('Dpp', data=dpp)
        
        
        
    