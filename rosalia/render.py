import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera
import matplotlib.animation as animation
import healpy as hp
from tqdm import tqdm
from IPython.display import Image

def render_healpix_movie(healpix_list, outname, figsize=(20,10), verbose=True):

    fig, ax = plt.subplots(1, figsize=figsize)
    camera = Camera(fig)

    for i in tqdm(range(len(healpix_list))):
        plt.axes(ax)
        hp.cartview(healpix_list[i], nest=True, min=0.07, max=0.2, hold=True)
        camera.snap()

    if verbose: print("Saving movie...")
    movie = camera.animate()
    writer = animation.PillowWriter(fps=30,
                                 metadata=dict(artist='Me'))

    movie.save(outname, writer=writer)

    if verbose: Image(outname)
    return(outname)
