"""The seven-dilate/six-fill FIJI LMN sequence and central component choice."""

import math

import numpy as np
import scipy.ndimage as ndi


def prepare_lmn(labels, dilation_iterations=7, hole_fill_iterations=6):
    if (dilation_iterations, hole_fill_iterations) != (7, 6):
        raise ValueError('Only the validated seven-dilate/six-fill sequence is supported')
    mask = labels > 0
    square = np.ones((1,3,3), bool)
    cross = ndi.generate_binary_structure(2, 1)
    for step in range(7):
        mask = ndi.binary_dilation(mask, structure=square)
        if step >= 1:
            mask = np.stack([ndi.binary_fill_holes(plane, structure=cross)
                             for plane in mask])
    return mask


def select_central_lmn(prepared, settings):
    components, number = ndi.label(prepared, structure=np.ones((3,3,3), bool))
    if number == 0:
        raise ValueError('No plausible LMN nucleus survived segmentation')
    counts = np.bincount(components.ravel())
    centers = ndi.center_of_mass(np.ones(prepared.shape, np.uint8), components,
                                 index=list(range(1, number+1)))
    xy_center = (np.array(prepared.shape[1:])-1)/2
    candidates = []
    for component, center in enumerate(centers, 1):
        if counts[component]:
            distance = float(np.linalg.norm(np.asarray(center[1:])-xy_center))
            candidates.append((distance, -int(counts[component]), component))
    candidates.sort()
    best = candidates[0]
    gap = candidates[1][0]-best[0] if len(candidates)>1 else math.inf
    confidence = ('high' if gap >= float(settings['center_gap_high_px']) else
                  'medium' if gap >= float(settings['center_gap_medium_px']) else 'low')
    return components == best[2], dict(selected_lmn_component=best[2],
        lmn_components_detected=number, selection_confidence=confidence,
        selected_lmn_volume_voxels=int(counts[best[2]]),
        center_distance_px=best[0], next_center_gap_px=gap)
