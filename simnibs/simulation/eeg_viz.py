from typing import Union

import numpy as np
import pyvista as pv

from simnibs.simulation.eeg import HEMISPHERES, FsAverage


class FsAveragePlotter:
    def __init__(self, subdivision: Union[int, str] = 7, surface: str = "inflated"):
        self.fsavg = FsAverage(subdivision)
        self.brain = self.surface_as_multiblock(surface)
        if surface == "inflated":
            self.brain["lh"].points[:, 0] -= np.abs(self.brain["lh"].points[:, 0].max())
            self.brain["rh"].points[:, 0] += np.abs(self.brain["rh"].points[:, 0].min())
        self.overlays = self.brain.copy()

    def surface_as_multiblock(self, surface):
        """Return the specified surface as a PyVista MultiBlock object."""
        surf = (
            self.fsavg.get_central_surface()
            if surface == "central"
            else self.fsavg.get_surface(surface)
        )
        mb = pv.MultiBlock()
        for hemi in surf:
            mb[hemi] = pv.make_tri_mesh(surf[hemi]["points"], surf[hemi]["tris"])
        return mb

    def add_curvature(self):
        curv = self.fsavg.get_morph_data("curv")
        for h in curv:
            self.brain[h]["curv"] = curv[h].astype(float)
            self.brain[h]["curv_bin"] = np.where(curv[h] > 0, 1 / 3, 2 / 3)

    def add_overlay(self, data, name: str):
        for h in HEMISPHERES:
            self.overlays[h][name] = data[h].astype(float)

    def remove_overlay(self, name: str):
        for h in HEMISPHERES:
            self.overlays[h].point_data.remove(name)

    def remove_all_overlays(self):
        for h in HEMISPHERES:
            self.overlays[h].clear_data()

    def set_active_overlay(self, name: str):
        assert name is not None
        for h in HEMISPHERES:
            self.overlays[h].set_active_scalars(name)

    def apply_threshold(self, threshold, use_abs: bool = False):
        overlay_thres = pv.MultiBlock()
        for h in HEMISPHERES:
            point_data = self.overlays[h].active_scalars
            cell_data = point_data[self.brain[h].faces.reshape(-1, 4)[:, 1:]].mean(-1)
            cell_data = np.abs(cell_data) if use_abs else cell_data
            overlay_thres[h] = self.overlays[h].remove_cells(cell_data < threshold)
        return overlay_thres

    def plot(
        self,
        overlay=None,
        threshold=None,
        use_abs=False,
        name=None,
        brain_kwargs=None,
        overlay_kwargs=None,
        plotter_kwargs=None,
        plotter=None,
    ):
        name = name or "temporary scalars"
        kw_brain = {}  # dict(
        #            scalars="curv_bin", cmap="gray", clim=(0, 1), show_scalar_bar=False
        #        )
        if brain_kwargs:
            kw_brain.update(brain_kwargs)
        kw_overlay = dict(
            annotations={threshold: "Threshold"} if threshold else None,
            # cmap="jet" if use_abs else "hot",
            # clim=None if use_abs else (0, None),
        )
        if overlay_kwargs:
            kw_overlay.update(overlay_kwargs)
        kw_plotter = plotter_kwargs or {}

        p = plotter or pv.Plotter(**kw_plotter)
        # Only plot background mesh if there is a "texture" or the overlay is
        # thresholded or there is no overlay
        if (
            all(h.active_scalars is not None for h in self.brain)
            or threshold
            or overlay is None
        ):
            p.add_mesh(self.brain, **kw_brain)

        if overlay is not None:
            if isinstance(overlay, str):
                self.set_active_overlay(overlay)
                remove_overlay = False
            else:
                # Temporarily add as overlay
                self.add_overlay(overlay, name)
                self.set_active_overlay(name)
                remove_overlay = True

            over = (
                self.apply_threshold(threshold, use_abs)
                if threshold
                else self.overlays.copy()
            )
            p.add_mesh(over, **kw_overlay)

            if remove_overlay:
                self.remove_overlay(name)

        p.view_xy(True)
        # p.camera.zoom(np.sqrt(2))
        p.camera.zoom(1.3)

        return p
