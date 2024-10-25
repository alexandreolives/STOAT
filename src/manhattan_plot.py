from itertools import cycle
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt

def manhattanplot(data, chrom="CHR", pos="POS", pv="P", snp="ID", 
                  logp=True, ax=None,
                  marker=".", color="#3B5488,#53BBD5", alpha=0.8,
                  title=None, xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$",
                  xtick_label_set=None, CHR=None, xticklabel_kws=None,
                  suggestiveline=1e-5, genomewideline=5e-8, sign_line_cols="#D62728,#2CA02C", hline_kws=None,
                  sign_marker_p=None, sign_marker_color="r",
                  is_annotate_topsnp=False, highlight_other_SNPs_indcs=None,
                  highlight_other_SNPs_color="r", highlight_other_SNPs_kwargs=None,
                  text_kws=None, ld_block_size=50000, **kwargs):

    data[[chrom]] = data[[chrom]].astype(str)  # make sure all the chromosome id are character.

    # Draw the plot and return the Axes
    if ax is None:
        # ax = plt.gca()
        _, ax = plt.subplots(figsize=(9, 3), facecolor="w", edgecolor="k")  # default

    if xticklabel_kws is None:
        xticklabel_kws = {}
    if hline_kws is None:
        hline_kws = {}
    if text_kws is None:
        text_kws = {}

    if "," in color:
        color = color.split(",")
    colors = cycle(color)

    last_xpos = 0
    xs_by_id = []  # use for collecting chromosome's position on x-axis
    x, y, c = [], [], []
    sign_snp_sites = []
    for seqid, group_data in data.groupby(by=chrom, sort=False):  # keep the raw order of chromosome

        if (CHR is not None) and (seqid != CHR):
            continue

        color = next(colors)
        for i, (site, p_value) in enumerate(zip(group_data[pos], group_data[pv])):
            if p_value == 0:
                p_value = 1e-300  # set it to a very small value if p-value is 0.

            y_value = -np.log10(p_value) if logp else p_value
            x.append(last_xpos + site)
            y.append(y_value)
            c.append(sign_marker_color if ((sign_marker_p is not None) and (p_value <= sign_marker_p)) else color)

            if (snp is not None) and (sign_marker_p is not None) and (p_value <= sign_marker_p):
                snp_id = group_data[snp].iloc[i]
                sign_snp_sites.append([last_xpos + site, y_value, snp_id])  # x_pos, y_value, text

        xs_by_id.append([seqid, last_xpos + (group_data[pos].iloc[0] + group_data[pos].iloc[-1]) / 2])
        last_xpos = x[-1]  # keep track so that chromosome will not overlap in the plot.

    if not x:
        raise ValueError("zero-size array to reduction operation minimum which has no "
                         "identity. This could be caused by zero-size array of ``x`` "
                         "in the ``manhattanplot(...)`` function.")

    if "marker" not in kwargs:
        kwargs["marker"] = marker

    # plot the main manhattan dot plot
    ax.scatter(x, y, c=c, alpha=alpha, edgecolors="none", **kwargs)

    highlight_other_SNPs_kwargs = dict() if highlight_other_SNPs_kwargs is \
                                            None else highlight_other_SNPs_kwargs

    # highlight other SNPs
    if highlight_other_SNPs_indcs is not None:
        for i in highlight_other_SNPs_indcs:
            ax.scatter(x[i], y[i], c=highlight_other_SNPs_color,
                       alpha=alpha, edgecolors="none", **highlight_other_SNPs_kwargs)

    # Add GWAS significant lines
    if "color" in hline_kws:
        hline_kws.pop("color")

    sign_line_cols = sign_line_cols.split(",") if "," in sign_line_cols else sign_line_cols
    if suggestiveline is not None:
        ax.axhline(y=-np.log10(suggestiveline) if logp else suggestiveline, color=sign_line_cols[0], **hline_kws)
    if genomewideline is not None:
        ax.axhline(y=-np.log10(genomewideline) if logp else genomewideline, color=sign_line_cols[1], **hline_kws)

    if CHR is None:

        if xtick_label_set is not None:
            ax.set_xticks([v for c, v in xs_by_id if c in xtick_label_set])
            ax.set_xticklabels([c for c, v in xs_by_id if c in xtick_label_set], **xticklabel_kws)
        else:
            ax.set_xticks([v for c, v in xs_by_id])
            ax.set_xticklabels([c for c, v in xs_by_id], **xticklabel_kws)

    else:
        ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.set_xlim(0, x[-1])
    ax.set_ylim(ymin=min(y), ymax=1.2 * max(y))

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return ax
