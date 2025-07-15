import numpy as np

import xtrack as xt

def patch_aperture_with(missing, patch, missing_apertures, line):
    if isinstance(missing, str) or not hasattr(missing, '__iter__'):
        missing = [missing]
    for nn in missing:
        if nn not in missing_apertures:
            print(f"No need to patch {nn}, aperture is present")
            continue
        if nn not in line.element_names:
            print(f"Element {nn} not found in line! Skipping aperture patching..")
            continue
        if isinstance(patch, str):
            if patch not in line.element_names:
                raise ValueError("Could not find patch aperture!")
            patch = line[patch].copy()
        line.insert_element(index=nn, element=patch,
                            name=nn+'_aper_patch')