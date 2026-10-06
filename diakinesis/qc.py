"""Compact labeled QC overlay: unique retained colors, red excluded, white LMN."""

import colorsys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi


COLORS = ['#2474b5','#25a542','#e2bc19','#16b7c7','#bd52c9',
          '#e58a2d','#6859b5','#788922','#4b9b76','#c48635',
          '#565dca','#27a9a1','#9e6cbb','#adb32d','#4c87d9']


def distinct_colors(count):
    """Stable within-stack palette, extending beyond the common 4–8 objects."""
    if count <= len(COLORS):
        return COLORS[:count]
    extra=[]
    for i in range(count-len(COLORS)):
        hue=.08+.74*((i*.618033988749895)%1)
        saturation=.62 if i%2 else .88
        value=.82 if (i//2)%2 else .98
        extra.append(matplotlib.colors.to_hex(colorsys.hsv_to_rgb(hue,saturation,value)))
    return COLORS+extra


def save_qc(path, source, dapi_image, labels, lmn_mask, objects, bins, confidence):
    projection = np.max(dapi_image,axis=0)
    lo,hi = np.percentile(projection,[1,99.5])
    background=np.clip((projection-lo)/max(hi-lo,1e-12),0,1)
    rgb=np.repeat((background*.65)[:,:,None],3,axis=2)
    selected=[o for o in objects if o['selected']]
    excluded=[o for o in objects if not o['selected']]
    palette=distinct_colors(len(selected))
    color_for={o['dapi_object_id']:palette[i] for i,o in enumerate(selected)}
    for obj in objects:
        oid=obj['dapi_object_id']
        mask=np.any(labels==oid,axis=0)
        color=matplotlib.colors.to_rgb(color_for[oid] if obj['selected'] else '#e31a1c')
        rgb[mask]=.25*rgb[mask]+.75*np.asarray(color)
    projection_lmn=np.any(lmn_mask,axis=0)
    edge=projection_lmn & ~ndi.binary_erosion(projection_lmn)
    rgb[edge]=(1,1,1)
    fig,ax=plt.subplots(figsize=(5.6,5.8),constrained_layout=True)
    ax.imshow(rgb,origin='lower')
    for obj in objects:
        oid=obj['dapi_object_id']
        color=color_for[oid] if obj['selected'] else '#e31a1c'
        ax.text(obj['centroid_x_px'],obj['centroid_y_px'],str(oid),
                color='white',ha='center',va='center',fontsize=6.5,weight='bold',
                bbox=dict(boxstyle='round,pad=.12',facecolor=color,edgecolor='black',alpha=.9))
    ax.set_title(f'{source}\nDAPI bins={bins}   Retained={len(selected)}   '
                 f'Excluded={len(excluded)}   LMN confidence={confidence}',fontsize=9)
    ax.set_xticks([]);ax.set_yticks([])
    fig.savefig(path,dpi=120)
    plt.close(fig)
