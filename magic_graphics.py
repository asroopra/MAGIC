#!/usr/bin/env python
"""
Drawing module
"""
import numpy as np

from plotly.subplots import make_subplots
import plotly.express as px
from pathlib import Path

import pandas as pd

########################################################################################################


def histplot(summary_df, chip_vals, padj_cutoff, fig_folder):
    """
    chip_vals = dict of dfs key = expt, val = [master df, query df]
    where master and query dfs are genes and chip vals as columns
    """

    # get list of experiments to plot (padj<cutoff)
    expts = summary_df[summary_df["padj"] <= padj_cutoff]["Experiment"].to_list()

    # cycle through experiments

    for ex in expts:
        # get TF name
        try:
            tf = ex.split(":")[1]
        except:
            tf = ex.split(":")[0]

        # get master and query dfs.
        mdf = chip_vals[ex][0]
        qdf = chip_vals[ex][1]

        # combine dfs (place one on top of other)
        df = pd.concat([mdf, qdf])

        # Add 'type' column
        df["type"] = ["master"] * len(mdf) + ["query"] * len(qdf)

        # get gene names for annotation
        genes = df["GENE"].to_list()

        # D argument for plotting
        argD = float(summary_df[summary_df["Experiment"] == ex]["argD"])

        # make figure with 3 subplots - in a column with shared x axis
        fig = make_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            row_heights=[0.4, 0.1, 0.5],
            vertical_spacing=0.02,
        )

        # make CDFs from df.
        """
		NOTE: px.ecdf and px.histogram return lists of dicts.
		For histograms with marginal plot:
		  the hist data is accessed in slot [0]     eg fig.data[0]
		  the marginal data is accessed in slot [1].eg fig.data[1]
		If there are multiple hists as here (2 hists):
		  the hist data is accessed in slot [0] and [2]
		  the marginal data is accessed in slot [1] and [3].	
		so the hist and marg data is stored as pairs		
		"""

        fig_ecdf = px.ecdf(df, x=ex, color="type")

        fig_hist = px.histogram(
            df,
            x=ex,
            color="type",
            histnorm="probability density",
            hover_name=genes,
            opacity=0.5,
            marginal="rug",
        )

        # add cdf to top panel
        fig.add_trace(fig_ecdf.data[0], row=1, col=1)
        fig.add_trace(fig_ecdf.data[1], row=1, col=1)
        fig.add_vline(argD)

        # add rugs to middle panel
        fig.add_trace(fig_hist.data[1], row=2, col=1)
        fig.add_trace(fig_hist.data[3], row=2, col=1)
        fig.add_vline(argD)

        # add pdfs to bottom panel
        fig.add_trace(fig_hist.data[0], row=3, col=1)
        fig.add_trace(fig_hist.data[2], row=3, col=1)
        fig.add_vline(argD)

        # place x axis title on bottom of bottom subpplot
        fig.update_xaxes(title_text="ChIP signal", row=3, col=1)

        # y axis titles
        fig.update_yaxes(title_text="Cumulative", row=1, col=1)
        fig.update_yaxes(title_text="Probability", row=3, col=1)

        fig.update_layout(barmode="overlay")

        # make title.  get padj to 4 decimal.
        padj = float(summary_df[summary_df["Experiment"] == ex]["padj"])
        title = f"Factor:{tf}  padj={padj:.2e}"  # {x:.2e} return 2 decimal places to scientific notation numbers

        fig.update_layout(title_text=title)

        f_path = str(Path(f"{fig_folder}/{tf}.html"))

        fig.write_html(f_path)


########################################################################################################
def summary_fig(summary_df, padj_cutoff, fig_path):
    # sort by df padj and trim to significant factors
    if len(summary_df) > 1:
        summary_df = summary_df.sort_values(by="Score", ascending=True)

    df = summary_df[summary_df["padj"] <= padj_cutoff].copy()

    fig = px.bar(
        df, x="Score", y="TF", color="padj", color_continuous_scale="bluered_r"
    )
    fig.write_html(fig_path)
