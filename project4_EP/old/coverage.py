import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def load_maps(desi_fp, act_fp):
    desi = hp.read_map(desi_fp)
    act = hp.read_map(act_fp)
    return desi + 2*act


def plotter(maps):
    cmap = ListedColormap(["whitesmoke", "yellow", "steelblue",
                           "mediumseagreen"], "overlap")
    hp.mollview(maps, coord=['C', 'G'], fig=2, cmap=cmap, max=3, cbar=False,
                notext=True, title=None)
    hp.graticule(alpha=0.5)
    hp.projtext(-0.9, -0.5, 'DESI', size=20)
    hp.projtext(8.9, -10., 'Intersection', size=20)
    hp.projtext(9.1, 12, 'advACT', size=20)
    plt.savefig('coverages.png', dpi=1080)


if __name__ == '__main__':
    maps = load_maps('project4_EP/overlaps/masks/DESI.fits',
                     'project4_EP/overlaps/masks/AdvACT.fits')
    plotter(maps)
