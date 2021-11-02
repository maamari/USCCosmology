
import os
import streamlit as st
import numpy as np
import healpy as hp
from matplotlib import cm
from matplotlib.colors import ListedColormap
import pandas as pd

st.set_option('deprecation.showPyplotGlobalUse', False)

cmap = ListedColormap(["tab:grey", "tab:blue", "tab:red", "tab:purple"], "overlap")

cmb_list = sorted([
            "Planck-Gal-70",
            "SPTPol",
            "SPTSZ",
            "ACTPol",
            "AdvACT",
            "Simons Observatory (Preliminary)",
            "CMB-S4 (Preliminary)",
            ])

galaxy_list = sorted([
               "BOSS-DR10",
               "BOSS-North",
               "DES",
               "DESI",
               "LSST",
               "CMASS-North",
               "eBOSS-North",
               "LOWZ-North",
               ])

cluster_catalogs = sorted([
                    "redMaPPer",
                    "WHL",
                    "DR9",
                    "PSZ2",
                    "AdvACT",
                    ])

def get_fsky(input_mask, threshold=0.1):
    """get the fraction of the observable sky

    Parameters
    ---------
    input_mask: np.ndarray
        healpy array indicating the input mask (0: masked, 1: visible)
    threshold: int
        mask cutoff value
    """
    return len(input_mask[input_mask > threshold]) / len(input_mask)

if __name__ == "__main__":
    st.title("CMB x Optical/IR Survey Overlaps")
    st.markdown("*Visualize overlaps between Cosmic Microwave Background (CMB) and Optical/IR Surveys*")

    # --------------
    # CMB Experiment
    # --------------

    # add a checkbox to select the CMB experiment
    st.sidebar.markdown("# CMB")
    cmb = st.sidebar.selectbox("(blue)", cmb_list, index=cmb_list.index("Planck-Gal-70"))
    cmb_fname = os.path.join(".", "masks", cmb + ".fits")
    cmb_mask = hp.read_map(cmb_fname)

    # add button for adding foregrounds (optional)
    add_foregrounds = st.sidebar.checkbox("Add Planck Foregrounds")
    foregrounds_fname = os.path.join(".", "masks", "Planck-Gal-70.fits")

    if add_foregrounds:
        foregrounds = hp.read_map(foregrounds_fname)
        cmb_mask = np.minimum(cmb_mask, foregrounds)

    # calculate and print the CMB f_sky
    cmb_fsky = get_fsky(cmb_mask)
    st.sidebar.text(f"f_sky = {cmb_fsky:.2f}")

    st.sidebar.markdown("---")

    # -------------
    # Galaxy Survey
    # -------------

    # add a checkbox to select the galaxy survey
    st.sidebar.markdown("# Optical/IR")
    galaxy = st.sidebar.selectbox("(red)", galaxy_list, index=galaxy_list.index("DESI"))
    gal_fname = os.path.join(".", "masks", galaxy+".fits")
    gal_mask = hp.read_map(gal_fname)

    # calculate and print the galaxy f_sky
    gal_fsky = get_fsky(gal_mask)
    st.sidebar.text(f"f_sky = {gal_fsky:.2f}")

    st.sidebar.markdown("---")

    # -------------
    # Cluster catalog
    # -------------

    # add a checkbox to select the galaxy survey
    st.sidebar.markdown("# Cluster catalog")
    catalog = st.sidebar.selectbox("(grey)", cluster_catalogs, index=cluster_catalogs.index("DR9"))
    catalog_fname = os.path.join(".", "masks", catalog+".csv")
    catalogs_in = pd.read_csv(catalog_fname,header=None)
    catalogs_in.columns = ['glon','glat']

    # calculate and print the galaxy f_sky
    num_clusters = len(catalogs_in['glon'])
    st.sidebar.text(f"Number of clusters = {num_clusters}")

    # plot the CMB (pixel values=1) and galaxy (pixel values=2) masks
    hp.mollview(cmb_mask,coord=['C','G'], fig=1, sub=121, cmap=cmap, title=cmb, cbar=False, max=3,notext=True)
    hp.mollview(2 * gal_mask,coord=['C','G'], fig=1, sub=122, cmap=cmap, title=galaxy, cbar=False, max=3, notext=True)
    st.pyplot()

    # calculate the overlap
    overlap_mask = cmb_mask + 2 * gal_mask
    overlap_only_mask = 3 * np.logical_and(cmb_mask, gal_mask)

    if st.checkbox("Show overlap only"):
        overlap_mask = overlap_only_mask

    overlap_fsky = get_fsky(overlap_only_mask)
    st.text(f"Overlap f_sky = {overlap_fsky:.2f}")

    hp.mollview(overlap_mask,coord=['C','G'], fig=2, cmap=cmap, title=f"{cmb} x {galaxy}", cbar=False, max=3)

    if st.checkbox(f"Show {catalog} clusters"):
        hp.visufunc.projscatter(catalogs_in['glon'],catalogs_in['glat'],c='black',s=0.01,lonlat= True)
    st.text(f"Number of clusters = {num_clusters}")
    st.pyplot()

    st.sidebar.markdown("---")
