import numpy as np
import matplotlib.pyplot as plt
import configparser
import sys
import h5py
from matplotlib.lines import Line2D

gPI = 3.141592653589793238462
gD2R = gPI / 180.0
gC = 1
gE0 = 0.511875
gME = gE0 / (gC * gC)
gRE = 6371000
T0 = 1.3802
T1 = 0.7405
E0=0.511875
c = 1.0
me = E0 / c**2
BE = 0.312 # T

K1 = 0.35
mu1 = 350

K2 = 0.45
mu2 = 500


def fname_base(run_id):
    return '../output/' + run_id + '/' + run_id


def parse_ini(run_id):
    config = configparser.ConfigParser()
    ini_fname= fname_base(run_id) + '.ini'

    config.read(ini_fname)

    paras = {}

    try:
        paras['E_MIN'] = float(config.get('basic', 'Emin'))  # MeV
        paras['E_MAX'] = float(config.get('basic', 'Emax'))  # MeV
        paras['nalpha0'] = int(config.get('basic', 'nalpha0'))
        paras['nE'] = int(config.get('basic', 'nE'))
        paras['nL'] = int(config.get('basic', 'nL'))
    except Exception as e:
        print("section_name or option_name wrong, check the input file.")
        sys.exit(1)

    return paras


def Y(alpha0):
    y = np.sin(alpha0)
    return 2 * T0 * (1-y) + (T0 - T1) * (y * np.log(y) + 2 * y - 2 * np.sqrt(y))


def p2e(p, E0=0.511875, cv=1):
    return np.sqrt(p**2* cv**2 + E0**2) - E0


def interpAlpha(K, mu, L_):
    # for certain K and mu interpolation
    alphaRange = np.linspace(5*np.pi/180, np.pi/2, 1000)
    L = (Y(alphaRange) * np.sqrt(BE) / (np.sin(alphaRange) * K))**2
    L = L[::-1]
    
    alpha0 = np.interp(L_, L, alphaRange)
    alpha0 = 95 * np.pi / 180 - alpha0
    # alpha0 = alpha0[::-1] #ERROR here
    p = np.sqrt(mu * BE * 2 * me / (L_**3 * np.sin(alpha0)**2))
    logE = np.log(p2e(p))
    return alpha0, logE


def alpha2i(alpha0):
    dalpha0 = (np.pi / 2.0) / paras['nalpha0']
    return (alpha0 - dalpha0 / 2.0) / dalpha0
    

def logE2j(logE):
    dlogE = (np.log(paras['E_MAX']) - np.log(paras['E_MIN'])) / paras['nE']
    return (logE - np.log(paras['E_MIN']) - dlogE / 2.0) / dlogE


if __name__ == '__main__':
    
    run_id = 'Radial'
    paras = parse_ini(run_id)
    
    hdL = (8 - 3) / (2 * paras['nL'])
    L_ = np.linspace(3+hdL, 8-hdL, paras['nL'])
    y_steady = 8 ** 7 / (8 ** 7 - 3 ** 7) - (3 ** 7 * 8 ** 7) / ((8 ** 7 - 3 ** 7) * L_ ** 7)

    a1, logE1 = interpAlpha(K1, mu1, L_)
    a2, logE2 = interpAlpha(K2, mu2, L_)

    file = h5py.File(fname_base(run_id) + '_data.h5', 'r')
    
    f01 = file['f/1'][:,:,:]
    f10 = file['f/10'][:,:,:]
    f500 = file['f/500'][:, :, :]
    
    f01_1 = []
    f10_1 = []
    f500_1 = []
    
    for k in range(paras['nL']):
        rawi = alpha2i(a1[k])
        rawj = logE2j(logE1[k])
        i = int(np.floor(rawi))
        j = int(np.floor(rawj))
        wi = rawi - i
        wj = rawj - j
        value01 = f01[i][j][k] * (1-wi) * (1-wj) + f01[i+1][j][k] * wi * (1-wj) + f01[i][j+1][k] * (1-wi) * wj + f01[i+1][j+1][k] * wi * wj
        value10 = f10[i][j][k] * (1-wi) * (1-wj) + f10[i+1][j][k] * wi * (1-wj) + f10[i][j+1][k] * (1-wi) * wj + f10[i+1][j+1][k] * wi * wj
        value500 = f500[i][j][k] * (1 - wi) * (1 - wj) + f500[i + 1][j][k] * wi * (1 - wj) + f500[i][j + 1][k] * (
                    1 - wi) * wj + f500[i + 1][j + 1][k] * wi * wj
        f01_1.append(value01)
        f10_1.append(value10)
        f500_1.append(value500)

    f01_2 = []
    f10_2 = []
    f500_2 = []

    for k in range(paras['nL']):
        rawi = alpha2i(a2[k])
        rawj = logE2j(logE2[k])
        i = int(np.floor(rawi))
        j = int(np.floor(rawj))
        wi = rawi - i
        wj = rawj - j
        value01 = f01[i][j][k] * (1-wi) * (1-wj) + f01[i+1][j][k] * wi * (1-wj) + f01[i][j+1][k] * (1-wi) * wj + f01[i+1][j+1][k] * wi * wj
        value10 = f10[i][j][k] * (1-wi) * (1-wj) + f10[i+1][j][k] * wi * (1-wj) + f10[i][j+1][k] * (1-wi) * wj + f10[i+1][j+1][k] * wi * wj
        value500 = f500[i][j][k] * (1 - wi) * (1 - wj) + f500[i + 1][j][k] * wi * (1 - wj) + f500[i][j + 1][k] * (
                    1 - wi) * wj + f500[i + 1][j + 1][k] * wi * wj
        f01_2.append(value01)
        f10_2.append(value10)
        f500_2.append(value500)

    file1d = h5py.File('../output/Radial1d/FVM1D.h5', 'r')
    
    L = file1d['L'][:]
    fmat = file1d['f'][:, :]

    fig = plt.figure(figsize=(20, 8))

    ax1 = fig.add_subplot(1, 2, 1)

    l1_1, = plt.plot(L, fmat[0], color='black', linestyle=':', label="initial condition")
    l1_2, = plt.plot(L, fmat[1], 'x', color='C0', label="1D FVM", markevery=2,  markeredgewidth=2)
    l2_2, = plt.plot(L_, f01_1, color='C0', label="Sayram")
    l1_3, = plt.plot(L, fmat[10], 'x',  color='C1', markevery=2, markeredgewidth=2)
    l2_2, = plt.plot(L_, f10_1, color='C1')
    l3, = plt.plot(L, y_steady, 'o', color='black', label="steady state", mfc='white', ms=10, markevery=2,  markeredgewidth=2)
    l3_1, = plt.plot(L_, f500_1, color='C2')
    l3_2, = plt.plot(L_, fmat[-1], 'x', color='C2', markevery=2,  markeredgewidth=2)

    proxy_lines1 = [Line2D([0], [0], linestyle=':', color='black'),
                    Line2D([0], [0], linestyle='none', marker='o', color='black', mfc='white', ms=10, markevery=2,  markeredgewidth=2),
                    Line2D([0], [0], linestyle='none', marker='x', color='black', markevery=2,  markeredgewidth=2),
                    Line2D([0], [0], color='black')]

    plt.legend(proxy_lines1, ["initial condition", "steady state", "1D FVM", "Sayram"])

    ax1.set_title(r'$\mu=$' + str(mu1) + ' MeV/G, K=' + str(K1) + ' G' + r'$^{1/2}$R' + r'$\mathrm{_E}$')

    plt.xlabel('L')
    plt.ylabel('f')

    ax2 = fig.add_subplot(1, 2, 2)

    l1_1 = plt.plot(L, fmat[0], color='black', linestyle=':')
    l1_2 = plt.plot(L, fmat[1], 'x', color='C0', markevery=2,  markeredgewidth=2)
    l2_2 = plt.plot(L_, f01_2, color='C0', label="    t = 100")
    l1_3 = plt.plot(L, fmat[10], 'x', color='C1', markevery=2,  markeredgewidth=2)
    l2_2 = plt.plot(L_, f10_2, color='C1', label="  t = 1000")
    l3 = plt.plot(L, y_steady, 'o', color='black', mfc='white', ms=10, markevery=2,  markeredgewidth=2)
    l3_1 = plt.plot(L_, fmat[-1], 'x', color='C2', markevery=2, markeredgewidth=2)
    l3_2 = plt.plot(L_, f500_1, color='C2', label="t = 50000")

    legend = plt.legend(handlelength=0, loc='best')

    legend.get_texts()[0].set_color('C0')
    legend.get_texts()[1].set_color('C1')
    legend.get_texts()[2].set_color('C2')

    ax2.set_title(r'$\mu=$' + str(mu2) + ' MeV/G, K=' + str(K2) + ' G' + r'$^{1/2}$R' + r'$\mathrm{_E}$')

    plt.xlabel('L')
    plt.ylabel('f')

    plt.tight_layout()
    plt.show()
