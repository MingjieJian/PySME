# -*- coding: utf-8 -*-
"""
Provide Plotting utility for Jupyter Notebook using Plot.ly
Can also be used just for Plot.ly, which will then generated html files
"""
import json
import logging
from base64 import b64encode

import numpy as np
import plotly.graph_objs as go
import plotly.offline as py
from plotly.io import write_image
from scipy.constants import speed_of_light

from .plot_colors import PlotColors

try:
    import ipywidgets as widgets
    from IPython import get_ipython
    from IPython.display import display

    cfg = get_ipython()
    in_notebook = cfg is not None
except (AttributeError, ImportError, ModuleNotFoundError):
    in_notebook = False

FigureClass = go.FigureWidget if in_notebook else go.Figure

logger = logging.getLogger(__name__)
clight = speed_of_light * 1e-3
fmt = PlotColors()

try:
    import htmlmin
except ImportError:
    logger.info("Install htmlmin for minified html files")
    htmlmin = None

if in_notebook:
    py.init_notebook_mode()


class FitPlot:
    """Plot the sme solve fit, as iterations pass along"""

    def __init__(self, wave, spec):
        self.fig = FigureClass()
        self.fig.layout["xaxis"]["title"] = "Wavelength [Å]"
        self.fig.layout["yaxis"]["title"] = "Intensity"

        self.fig.add_scatter(x=wave, y=spec, name="Observation")

    def add_synth(self, wave, synth, iteration=0):
        """add a scatter plot to the plot"""
        self.fig.add_scatter(x=wave, y=synth, name=f"Iteration {iteration}")


class FinalPlot:
    """Big plot that covers everything"""

    def __init__(
        self,
        sme,
        segment=0,
        orig=None,
        labels=None,
    ):
        self.sme = sme
        self.wave = sme.wave
        self.spec = sme.spec
        self.mask = sme.mask
        self.smod = sme.synth
        self.orig = orig
        self.labels = labels if labels is not None else {}
        if sme.telluric is not None and self.smod is not None:
            for i in range(sme.nseg):
                if len(self.smod[i]) != 0:
                    self.smod[i] = self.smod[i] * sme.telluric[i]
        self.nsegments = len(self.wave)
        self.segment = segment
        self.wran = sme.wran
        self.lines = sme.linelist
        self.vrad = sme.vrad
        self.vrad = [v if v is not None else 0 for v in self.vrad]

        self.mask_type = "good"

        data, annotations = self.create_plot(self.segment)
        self.annotations = annotations

        # Add segment slider
        steps = []
        for i in range(self.nsegments):
            step = {
                "label": f"Segment {i}",
                "method": "update",
                "args": [
                    {"visible": [v == i for v in self.visible]},
                    {
                        "title": f"Segment {i}",
                        "annotations": annotations[i],
                        "xaxis": {
                            "range": list(self.wran[i]),
                            "linecolor": "black",
                            "linewidth": 2.4,
                            "mirror": True,
                            "ticks": "outside",
                            "tickcolor": "black",
                        },
                        "yaxis": {
                            "autorange": True,
                            "linecolor": "black",
                            "linewidth": 2.4,
                            "mirror": True,
                            "ticks": "outside",
                            "tickcolor": "black",
                        },
                    },
                ],
            }
            steps += [step]

        layout = {
            "dragmode": "select",
            "plot_bgcolor": "#fff",
            "selectdirection": "h",
            "title": f"Segment {self.segment}",
            "xaxis": {
                "title": "Wavelength [Å]",
                "linecolor": "black",
                "linewidth": 2.4,
                "mirror": True,
                "ticks": "outside",
                "tickcolor": "black",
            },
            "yaxis": {
                "title": "Intensity",
                "linecolor": "black",
                "linewidth": 2.4,
                "mirror": True,
                "ticks": "outside",
                "tickcolor": "black",
            },
            "annotations": annotations[self.segment],
            "sliders": [{"active": self.segment, "steps": steps}],
            "legend": {"traceorder": "reversed"},
        }
        if in_notebook:
            self.fig = FigureClass(data=data, layout=layout)

            # add selection callback
            self.fig.data[0].on_selection(self.selection_fn)
        else:
            data = [go.Scattergl(d) for d in data]
            self.fig = FigureClass(data=data, layout=layout)

        # Add button to save figure
        if in_notebook:
            self.button_save = widgets.Button(description="Save")
            self.button_save.on_click(self.save)

            # Add buttons for Mask selection
            self.button_mask = widgets.ToggleButtons(
                options=["Good", "Bad", "Continuum", "Line"],
                description="Mask",
            )
            self.button_mask.observe(self.on_toggle_click, "value")

            self.widget = widgets.VBox([self.button_mask, self.button_save, self.fig])
            display(self.widget)

    def save(self, _=None, filename="SME.html", **kwargs):
        """save plot to html file"""
        # Here we "hack" plotly to use base64 encoded data
        # since we have a lot of data to look at
        # this reduces the file size by about a factor 3
        # zipping it would provide another factor 2
        if filename.endswith(".html"):
            self.fig.layout.dragmode = "zoom"
            fig = self.fig.to_plotly_json()
            # Convert the data to base64
            for i in range(len(fig["data"])):
                for ax in ["x", "y"]:
                    data = np.asarray(fig["data"][i][ax]).astype("float32")
                    b64data = b64encode(data).decode("utf8")
                    fig["data"][i][ax] = {"data": b64data, "dtype": "float32"}

            data = json.dumps(fig["data"])
            layout = json.dumps(fig["layout"])
            # This is the entire html as generated by self.fig.to_html
            # except we replace the data and layout manually
            # and transform it from base64 to TypedArrays
            html = (
                r"""<html><head><meta charset="utf-8" /></head><body><div>"""
                r"""<script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: 'local'};</script>"""
                r"""<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>"""
                r"""<div id="56dedbe2-a786-4b0b-8232-b68dfd7b9eb5" class="plotly-graph-div" style="height:100%; width:100%;"></div>"""
                r"""<script type="text/javascript">"""
                r"""window.PLOTLYENV=window.PLOTLYENV || {};"""
                r"""if (document.getElementById("56dedbe2-a786-4b0b-8232-b68dfd7b9eb5")) {"""
                f"data = {data};"
                f"layout = {layout};"
                r"""const FromBase64 = function (str) {return new Uint8Array(atob(str).split('').map(function (c) { return c.charCodeAt(0); }));};"""
                r"""for (let i = 0; i < data.length; i++){"""
                r"""let data_x = data[i].x.data;"""
                r"""let data_y = data[i].y.data;"""
                r"""data[i].x = new Float32Array(FromBase64(data_x).buffer);"""
                r"""data[i].y = new Float32Array(FromBase64(data_y).buffer);"""
                r"""}"""
                r"""Plotly.newPlot("56dedbe2-a786-4b0b-8232-b68dfd7b9eb5", data, layout, {"responsive": true})};"""
                r"""</script></div></body></html>"""
            )
            # Finally try to minimize if possible
            if htmlmin is not None:
                html = htmlmin.minify(html, remove_empty_space=True)
            with open(filename, "w") as f:
                f.write(html)

        else:
            write_image(self.fig, filename)

    def show(self):
        self.fig.show()

    def shift_mask(self, x, mask):
        """shift the edges of the mask to the bottom of the plot,
        so that the mask creates a shape with straight edges"""
        for i in np.where(mask)[0]:
            try:
                if mask[i] == mask[i + 1]:
                    x[i] = x[i - 1]
                else:
                    x[i] = x[i + 1]
            except IndexError:
                pass

        return x

    def create_mask_points(self, x, y, mask, value):
        """Creates the points that define the outer edge of the mask region"""
        mask = mask != value
        x = np.copy(x)
        y = np.copy(y)
        y[mask] = 0
        x = self.shift_mask(x, mask)
        return x, y

    def create_plot(self, current_segment):
        """Generate the plot componentes (lines and masks) and line labels"""

        annotations = {}
        visible = []
        data = []
        line_mask_idx = {}
        cont_mask_idx = {}

        for seg in range(self.nsegments):

            k = len(visible)
            line_mask_idx[seg] = k
            cont_mask_idx[seg] = k + 1

            # The order of the plots is chosen by the z order, from low to high
            # Masks should be below the spectra (so they don't hide half of the line)
            # Synthetic on top of observation, because synthetic varies less than observation
            # Annoying I know, but plotly doesn't seem to have good controls for the z order
            # Or Legend order for that matter

            if (
                self.mask is not None
                and self.spec is not None
                and self.wave is not None
            ):
                # Line mask
                x, y = self.create_mask_points(
                    self.wave[seg], self.spec[seg], self.mask[seg], 1
                )

                data += [
                    dict(
                        x=x,
                        y=y,
                        fillcolor=fmt["LineMask"]["facecolor"],
                        fill="tozeroy",
                        mode="none",
                        name=self.labels.get("linemask", "Line Mask"),
                        hoverinfo="none",
                        legendgroup=2,
                        visible=current_segment == seg,
                    )
                ]
                visible += [seg]

                # Cont mask
                x, y = self.create_mask_points(
                    self.wave[seg], self.spec[seg], self.mask[seg], 2
                )

                data += [
                    dict(
                        x=x,
                        y=y,
                        fillcolor=fmt["ContMask"]["facecolor"],
                        fill="tozeroy",
                        mode="none",
                        name=self.labels.get("contmask", "Continuum Mask"),
                        hoverinfo="none",
                        legendgroup=2,
                        visible=current_segment == seg,
                    )
                ]
                visible += [seg]

            if self.spec is not None:
                # Observation
                data += [
                    dict(
                        x=self.wave[seg],
                        y=self.spec[seg],
                        line={"color": fmt["Obs"]["color"]},
                        name=self.labels.get("obs", "Observation"),
                        legendgroup=0,
                        visible=current_segment == seg,
                    )
                ]
                visible += [seg]

            # Synthetic, if available
            if self.smod is not None:
                data += [
                    dict(
                        x=self.wave[seg],
                        y=self.smod[seg],
                        name=self.labels.get("synth", "Synthetic"),
                        line={"color": fmt["Syn"]["color"]},
                        legendgroup=1,
                        visible=current_segment == seg,
                    )
                ]
                visible += [seg]
            if self.orig is not None:
                data += [
                    dict(
                        x=self.wave[seg],
                        y=self.orig[seg],
                        name=self.labels.get("orig", "Original"),
                        line={"color": fmt["Orig"]["color"], "dash": "dash"},
                        legendgroup=1,
                        visible=current_segment == seg,
                    )
                ]
                visible += [seg]

            # mark important lines
            if self.lines is not None and len(self.lines) != 0:
                seg_annotations = []
                xlimits = self.wave[seg][[0, -1]]
                xlimits = xlimits * (1 - self.vrad[seg] / clight)
                lines = (self.lines.wlcent > xlimits[0]) & (
                    self.lines.wlcent < xlimits[1]
                )
                lines = self.lines[lines]

                # Filter out closely packaged lines of the same species
                # Threshold for the distance between lines
                wlcent, labels = [], []
                threshold = np.diff(xlimits)[0] / 100
                species = np.unique(lines["species"])
                for sp in species:
                    sp_lines = lines["wlcent"][lines["species"] == sp]
                    diff = np.diff(sp_lines)
                    sp_mask = diff < threshold
                    if np.any(sp_mask):
                        idx = np.where(np.diff(sp_mask))[0]
                        idx = [0, *idx, len(sp_mask) + 1]
                        idx = np.unique(idx)
                        for i, j in zip(idx[:-1], idx[1:]):
                            sp_wmid = np.mean(sp_lines[i:j])
                            sp_label = f"{sp} +{j-i-1}"
                            wlcent += [sp_wmid]
                            labels += [sp_label]
                    else:
                        for sp_wmid in sp_lines:
                            wlcent += [sp_wmid]
                            labels += [sp]

                wlcent = np.array(wlcent)
                labels = np.array(labels)

                # Keep only the 100 stongest lines for performance
                # if "depth" in lines.columns:
                #     lines.sort("depth", ascending=False)
                #     lines = lines[:20]
                # else:
                #     idx = np.random.choice(len(lines), 20, replace=False)
                #     lines = lines[idx]

                x = wlcent * (1 + self.vrad[seg] / clight)
                if self.smod is not None and len(self.smod[seg]) != 0:
                    y = np.interp(x, self.wave[seg], self.smod[seg])
                else:
                    y = np.interp(x, self.wave[seg], self.spec[seg])

                if self.smod is not None and len(self.smod[seg]) > 0:
                    ytop = np.max(self.smod[seg])
                elif self.spec is not None and len(self.spec[seg]) > 0:
                    ytop = np.max(self.spec[seg])
                else:
                    ytop = 1

                for i, line in enumerate(labels):
                    seg_annotations += [
                        {
                            "x": x[i],
                            "y": y[i],
                            "xref": "x",
                            "yref": "y",
                            "text": f"{labels[i]}",
                            "hovertext": f"{wlcent[i]}",
                            "textangle": 90,
                            "opacity": 1,
                            "ax": 0,
                            "ay": 1.2 * ytop,
                            "ayref": "y",
                            "showarrow": True,
                            "arrowhead": 7,
                            "xanchor": "left",
                        }
                    ]
                annotations[seg] = seg_annotations

        self.visible = visible
        self.line_mask_idx = line_mask_idx
        self.cont_mask_idx = cont_mask_idx

        return data, annotations

    def selection_fn(self, trace, points, selector):
        """Callback for area selection, changes the mask depending on selected mode"""
        self.segment = self.fig.layout["sliders"][0].active
        seg = self.segment

        xrange = selector.xrange
        wave = self.wave[seg]
        mask = self.mask[seg]

        # Choose pixels and value depending on selected type
        if self.mask_type == "good":
            value = 1
            idx = (wave > xrange[0]) & (wave < xrange[1]) & (mask == 0)
        elif self.mask_type == "bad":
            value = 0
            idx = (wave > xrange[0]) & (wave < xrange[1])
        elif self.mask_type == "line":
            value = 1
            idx = (wave > xrange[0]) & (wave < xrange[1]) & (mask != 0)
            print(np.count_nonzero(idx))
        elif self.mask_type == "cont":
            value = 2
            idx = (wave > xrange[0]) & (wave < xrange[1]) & (mask == 1)
        else:
            return

        # Apply changes if any
        if np.count_nonzero(idx) != 0:
            self.mask[seg][idx] = value

            with self.fig.batch_update():
                # Update Line Mask
                m = self.line_mask_idx[seg]
                x, y = self.create_mask_points(
                    self.wave[seg], self.spec[seg], self.mask[seg], 1
                )
                self.fig.data[m].x = x
                self.fig.data[m].y = y

                # Update Cont Mask
                m = self.cont_mask_idx[seg]
                x, y = self.create_mask_points(
                    self.wave[seg], self.spec[seg], self.mask[seg], 2
                )
                self.fig.data[m].x = x
                self.fig.data[m].y = y

    def on_toggle_click(self, change):
        """Callback for mask mode selector buttons"""
        change = change["new"]
        if change == "Good":
            self.set_mask_good()
        elif change == "Bad":
            self.set_mask_bad()
        elif change == "Continuum":
            self.set_mask_continuum()
        elif change == "Line":
            self.set_mask_line()

    def set_mask_good(self, _=None):
        """Called by clicking the 'good' mask button"""
        self.set_mask_type("good")

    def set_mask_bad(self, _=None):
        """Called by clicking the 'bad' mask button"""
        self.set_mask_type("bad")

    def set_mask_line(self, _=None):
        """Called by clicking the 'line' mask button"""
        self.set_mask_type("line")

    def set_mask_continuum(self, _=None):
        """Called by clicking the 'continuum' mask button"""
        self.set_mask_type("cont")

    def set_mask_type(self, type):
        """Changes the mask change mode and chooses the current interactive tool"""
        self.mask_type = type
        self.fig.layout["dragmode"] = "select"

    def add(self, x, y, label=""):
        """adds a scatter plot to the image, and makes the necessary changes in the slider"""
        self.fig.add_scatter(x=x, y=y, name=label, legendgroup=10)
        self.visible += [-1]

        # Update Sliders
        steps = []
        for i in range(self.nsegments):
            step_visible = [(v == i) or (v == -1) for v in self.visible]
            step = {
                "label": f"Segment {i}",
                "method": "update",
                "args": [
                    {"visible": step_visible},
                    {
                        "title": f"Segment {i}",
                        "annotations": self.annotations[i],
                        "xaxis": {"range": list(self.wran[i])},
                        "yaxis": {"autorange": True},
                    },
                ],
            }
            steps += [step]

        self.fig.layout["sliders"][0]["steps"] = steps
