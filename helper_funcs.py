import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
import scipy.stats
from scipy.stats import pearsonr
from PIL import Image

### Data downloads ###

IMG_FOLDER = "/Users/iboyle/Files/TargetDiscovery/PELO/PELO-dev/FINAL Raw WB Images/"

FIGSHARE_FILE_MAP = {
    "CRISPRGeneEffect":"43346616",
    "OmicsCNGene":"43346913",
    "OmicsExpressionProteinCodingGenesTPMLogp1":"43347204",
    "OmicsDefaultModelProfiles":"43346961",
    "OmicsCNSegmentsProfile":"43346952",
    "OmicsSomaticMutationsMatrixHotspot":"43347528",
    "OmicsSomaticMutationsMatrixDamaging":"43347516",
    "Model":"43746708",
    "OmicsSignatures":"46500361",
}

def download_figshare_file(name, **kwargs):
    '''
    Download files from figshare and save them locally to data folder.
    '''
    html = "https://ndownloader.figshare.com/files/"
    file_num = FIGSHARE_FILE_MAP[name]
    print(f"Fetching {name} from {html}{file_num}")
    return pd.read_csv(html+file_num, **kwargs).to_csv(f"data/{name}.csv")

def download_other_files(name):
    '''
    Download files not in figshare and save them locally to data folder.
    '''
    if name == "ccle_proteomics":
        try:
            pd.read_excel(
                r"https://pmc.ncbi.nlm.nih.gov/articles/instance/7339254/bin/NIHMS1569324-supplement-TS2.xlsx", 
                sheet_name=1, storage_options={'User-Agent': 'Mozilla/5.0'}
            ).to_csv("data/ccle_proteomics.csv", index=None)
        except Exception as e:
            print("tried to read 'ccle_proteomics' file from PMC but failed. Download supplementary table 2 (NIHMS1569324-supplement-TS2.xlsx) from PMID: 31978347, and save as data/ccle_proteomics.csv")
    elif name == "HGNC_complete_set":
        pd.read_csv(
            r"https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2024-01-01.txt", 
            sep="\t"
        ).to_csv("data/HGNC_complete_set.csv", index=None)

def prep_folders():
    for folder in ["data", "outputs", "figures"]:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
def read_in_data(name, redownload=False, **kwargs):
    '''
    Read in local version of data, or fetch if absent.
    '''
    path = f"data/{name}.csv"
    if not os.path.exists("data"):
        prep_folders()
    if not os.path.exists(path) or redownload:
        if name in FIGSHARE_FILE_MAP.keys():
            download_figshare_file(name, **kwargs)
        else:
            download_other_files(name)
    else:
        print(f"Fetching local copy of {name}")
    return pd.read_csv(path, **kwargs)

def munge_ccle_proteomics(ccle_proteomics, model_meta):
    '''
    Take the raw CCLE proteomics data and convert to ModelID x protein matrix.
    '''
    #index by gene symbol only and drop metadata columns
    ccle_proteomics = ccle_proteomics.set_index("Gene_Symbol").T.drop([
        'Protein_Id', 'Description', 'Group_ID', 'Uniprot', 'Uniprot_Acc'
    ])
    #map from stripped cell line name name (everything before _) to ModelIDs
    ccle_proteomics.index = ccle_proteomics.index.str.split("_").str[0]
    ccle_proteomics = ccle_proteomics.merge(
        model_meta.reset_index().set_index("StrippedCellLineName")[["ModelID"]], 
        left_index=True, right_index=True, how="inner"
    ).set_index("ModelID")
    #drop any cell lines mapping to multiple ModelIDs, which are likely problematic
    dup_cls = ccle_proteomics.loc[ccle_proteomics.index.duplicated(keep=False)]
    ccle_proteomics = pd.concat([
        ccle_proteomics.loc[~ccle_proteomics.index.isin(dup_cls.index)], 
        dup_cls.groupby(dup_cls.index).mean()
    ])
    return ccle_proteomics

### Figure styling ###

DOUBLE_IN = 7.05
SINGLE_IN = 3.54
SPECIAL_IN = 4.72

COLORS = ["#008F85", "#F07167",  "#88527F",  "#2D314E", "#5D686F", "#DC3232"]
PELO_GUIDES = ["#003D39", "#008F85", "#00CCBE"] 
POS_CON = "#F07167"
NEG_CON = "#88527F"
ACCENT = "#5D686F"
BLACK = "#000000"

def set_params():
    from matplotlib import rcParams
    rcParams['axes.spines.right'] = False
    rcParams['axes.spines.top'] = False
    rcParams['savefig.dpi'] = 500
    rcParams['savefig.transparent'] = False
    rcParams['font.family'] = 'Arial'
    rcParams['figure.dpi'] = 200
    rcParams["savefig.facecolor"] = (1, 1, 1.0, 0.2)
    
    rcParams['pdf.fonttype']=42
    rcParams['ps.fonttype'] = 42
    
    rcParams['axes.titlesize'] = 6
    rcParams['axes.labelsize'] = 6
    rcParams['font.size'] = 5
    rcParams['legend.fontsize'] = 5

    sns.set_palette(COLORS) 
set_params()

SHORT_ONCO_NAMES = {
    'Ovary/Fallopian Tube':'Ovary', 
    'Bladder/Urinary Tract':'Bladder', 
    'Soft Tissue':'Soft', 
    'Esophagus/Stomach':'Stomach',
    'Peripheral Nervous System':'PNS', 
    'Biliary Tract':'Biliary',
    'Head and Neck':'Head/Neck', 
    'Ampulla of Vater':'Ampulla',
}

HALLMARK_ABBREV = {
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE":"UPR",
    "HALLMARK_P53_PATHWAY":"p53",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB":"TNFa via NFkB",
    "HALLMARK_G2M_CHECKPOINT":"G2M",
    "HALLMARK_MITOTIC_SPINDLE":"Mitotic Spindle",
    "HALLMARK_E2F_TARGETS":"E2F",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION":"OxPhos",
    "HALLMARK_MYC_TARGETS_V1":"MYC (v1)",
    "HALLMARK_MTORC1_SIGNALING":"MTORC1",
}

## Statistics ##

def compare_it(classes, data, test=scipy.stats.ttest_ind, nan_policy='omit', alternative="less", plot=True, **plot_kws):
    '''
    Run a two-class statistical test (ie. t-test).

    Parameters:
        classes - pandas boolean Series defining two classes
        data - pandas DataFrame to run the comparison on, with index matching classes
        test - statistical test function from scipy.stats
        nan_policy - input to test which determines nan handling
        alternative - input to test which determines if left-tailed (less), right-tailed (greater), or two sided (two-sided)
        plot - boolean indicating whether to plot the results

    Outputs:
        pandas DataFrame with n (# samples), p (p-value), Q (false discovery value), Effect Size (mean difference)
    '''
    #get number of non-null values in each comparison
    n = data.loc[data.index.isin(classes.index)].notnull().sum()
    #find pvals
    pvals = pd.Series(
        test(
            data.loc[data.index.isin(classes.index[classes == True])],
            data.loc[data.index.isin(classes.index[classes == False])],
            alternative=alternative, nan_policy=nan_policy
        ).pvalue, index=data.columns
    )
    #run fdr to get qvals
    fdr = pd.Series(scipy.stats.false_discovery_control(pvals.dropna()), index=pvals.dropna().index)
    #compute mean difference (effect size)
    effect_size = (
        data.loc[data.index.isin(classes.index[classes == True])].mean() -
        data.loc[data.index.isin(classes.index[classes == False])].mean()
    ).rename("Mean Difference")
    if plot: 
        quick_scatter(
            x=effect_size,
            y=-np.log10(fdr).rename("-log$_{10}$(q-value)"),
            **plot_kws
        )
    return pd.DataFrame({"n":n, "p":pvals, "q":fdr, "Effect Size":effect_size})

def compare_it_continuous(x, data, alternative="two-sided", plot=True, **plot_kws):
    '''
    Run a continuous statistical test, aka Pearson correlation.

    Parameters:
        x - pandas Series of values to compare to
        data - pandas DataFrame to run the comparison on, with index matching x
        alternative - input to test which determines if left-tailed (less), right-tailed (greater), or two sided (two-sided)
        plot - boolean indicating whether to plot the results

    Outputs:
        pandas DataFrame with n (# samples), p (p-value), Q (false discovery value), Effect Size (Pearson's R)
    '''
    #helper to filter NAs before running Pearson
    def filter_then_pearson(y, x):
        x, y = x.dropna().align(y.dropna(), join="inner")
        n = x.shape[0]
        if y.nunique() > 1:
            return [n, scipy.stats.pearsonr(y, x, alternative=alternative).pvalue]
        return [n, np.nan]
    #get number of non-null values in each comparison and pvals
    res = data.apply(lambda y: pd.Series(filter_then_pearson(y, x), index=["n", "p"])).T
    #find fdr q-value
    res["q"] = pd.Series(scipy.stats.false_discovery_control(res["p"].dropna()), index=res["p"].dropna().index)
    #compute effect size (correlation)
    res["Effect Size"] = data.corrwith(x)
    if plot: 
        quick_scatter(
            x=res["Effect Size"],
            y=-np.log10(res["q"]).rename("-log$_{10}$(q-value)"),
            **plot_kws
        )
    return res

def ttest_vs_negcons(df, var_col, neg_cons, to_test, groupby=None, alternative='less', vals="Relative Viability"):
    '''
    Compute a Student's t-test on a pandas dataframe.

    Parameters:
        df - pandas DataFrame
        var_col - column in df to treat as variable defining neg controls and test
        neg_cons - negative controls in var_col
        to_test - test values in var_col
        groupby - column defining groups to perform the test within
        alternative - which scipy.stats.ttest_ind alternative to test
        vals - column in df containing values to be tested (often Relative Viability)

    Outputs:
        table containing the p/q-values, effect size (mean difference), and n (# samples in neg+test)
    '''
    #get replicates for neg and test associated with each group
    if groupby:
        #get test/neg con values in each group
        groups = df.groupby(groupby).apply(
            lambda x: (
                x.loc[x[var_col].isin(to_test), vals].dropna().tolist(),
                x.loc[x[var_col].isin(neg_cons), vals].dropna().tolist() 
            )
        )
        #run the t-test per group
        results = groups.apply(lambda x: pd.Series(scipy.stats.ttest_ind(*x, alternative=alternative), index=["statistic", "p-val"]))
        #FDR correct all results
        results["q"] = scipy.stats.false_discovery_control(results["p-val"])
        #add effect size and #
        results["Effect Size"] = groups.apply(lambda x: np.mean(x[0]) - np.mean(x[1]))
        results["n"] = groups.apply(lambda x: len(x[0]) + len(x[1]))
        
    #get replicates for neg and test
    else:
        test = df.loc[df[var_col].isin(to_test), vals].tolist()
        neg = df.loc[df[var_col].isin(neg_cons), vals].tolist()
        results = pd.Series(scipy.stats.ttest_ind(test, neg, alternative=alternative), index=["statistic", "p-val"])
        results["Effect Size"] = np.mean(test) - np.mean(neg)
        results["n"] = len(test)+len(neg)
    return results

def to_scientific(val, prefix="$p$="):
    '''
    Format a float (val) scientifically.
    '''
    if val == 0:
        return prefix+"0.000"
    elif val < 1e-3:
        return prefix+"{:.2e}".format(val).replace("e", "x10$^{")+"}$"
    else:
        return prefix+"{:.3f}".format(val)
        
def pearson(x, y):
    '''
    Compute Pearson's R between x and y then return as pandas Series.
    '''
    corr = pd.concat([pd.Series(x), pd.Series(y)], axis=1, join="inner")
    res = pearsonr(corr.iloc[:,0], corr.iloc[:,1])
    return pd.Series({"R":res.statistic, "p-val":res.pvalue})

def pearson_str(x, y):
    '''
    Compute Pearson's R between x and y then return as a well-formatted string.
    '''
    corr = pd.concat([pd.Series(x), pd.Series(y)], axis=1, join="inner")
    res = pearsonr(corr.iloc[:,0], corr.iloc[:,1])
    return f"R={res.statistic.round(3)}, {to_scientific(res.pvalue)}"

def add_labels_to_ax(ax, y_vals, sig_annots, y_adj=0, line_adj=(0.25, -0.01), linewidth=0.5,  fontsize=5, **kws):
    '''
    Add labels with lines to axes.
    '''
    xticks = pd.Series(ax.get_xticks(), [t.get_text() for t in ax.get_xticklabels()])
    for x, y in y_vals.items():
        if x in xticks.index:
            ax.text(x=x, y=y+y_adj, s=sig_annots[x], ha="center", fontsize=fontsize, **kws)
            ax.hlines(
                y=y+line_adj[1], xmin=xticks[x]-line_adj[0], xmax=xticks[x]+line_adj[0], 
                color=ACCENT, linewidth=linewidth
            )
        elif not isinstance(x, str):
            ax.text(x=x, y=y+y_adj, s=sig_annots[x], ha="center", fontsize=fontsize, **kws)
            ax.hlines(
                y=y+line_adj[1], xmin=x-line_adj[0], xmax=x+line_adj[0], 
                color=ACCENT, linewidth=linewidth
            )

def add_bracket(ax, label, x_groups, y, dx=0.23, dy=0.05, linewidth=0.5):
    '''
    Add brackets to axes.
    '''
    xticks = pd.Series(ax.get_xticks(), [t.get_text() for t in ax.get_xticklabels()])
    yt = y
    for j, xs in enumerate(x_groups):
        y = yt if len(x_groups) == 1 else yt - dy 
        xticks_curr = xticks.loc[xs]
        n_brackets = len(xs)-1
        if n_brackets == 0:
            x1 = xticks_curr.iloc[0]
            xs = [x1+dx,x1+dx]
            ys = [y-dy, y]
            ax.plot(xs, ys, color="k", linewidth=linewidth, clip_on=False)
        else:
            for i in range(n_brackets):
                x1, x2 = xticks_curr.iloc[i], xticks_curr.iloc[i+1]
                xs = [x1+dx,x1+dx, x2+dx, x2+dx]
                ys = [y-dy, y, y, y-dy]
                ax.plot(xs, ys, color="k", linewidth=linewidth, clip_on=False)
        if j < len(x_groups)-1:
            x1, x2 = xticks.loc[x_groups[j]].mean(), xticks.loc[x_groups[j+1]].mean()
            xs = [x1+dx,x1+dx, x2+dx, x2+dx]
            ys = [yt-dy, yt, yt, yt-dy]
            ax.plot(xs, ys, color="k", linewidth=linewidth, clip_on=False)
            ax.text((x1+x2)/2+dx, yt, label, va="bottom", ha="center", clip_on=False)
    if len(x_groups) == 1:
        x1, x2 = xticks.loc[x_groups[0]].min(), xticks.loc[x_groups[0]].max()
        ax.text((x1+x2)/2+dx, yt, label, va="bottom", ha="center", clip_on=False)

## Plotting ##

def axline(ax, dir, loc, color=ACCENT, alpha=0.5, label=None):
    '''
    Add a horizontal or vertical styled line to an axis.
    '''
    if dir == "h":
        ax.axhline(loc, color=color, alpha=alpha, linestyle="--", zorder=-1, label=label)
    elif dir == "v":
        ax.axvline(loc, color=color, alpha=alpha, linestyle="--", zorder=-1, label=label)
    else:
        raise Exception("unknown dir for axline")       

def quick_scatter(
    x, y, hue=None, style_hue=False, labels=None, figsize=(5,5), diagonal=False, hline=None, vline=None, 
    title=None, legend_out=True, hide_legend=False, text_kws=dict(ha="right"), arrowprops=dict(arrowstyle='->', color=ACCENT), 
    return_ax=False, ax=None, use_lims=None, adj_kws=dict(), **kwargs
):
    '''
    Produce a nicely styled scatterplot.

    Parameters:
        x - named pandas Series for x-axis values
        y - named pandas Series for y-axis values
        hue - named pandas Series for colors
        style_hue - boolean indicating whether to style colors with different markers
        labels - points in plot to label
        figsize - size of figure to produce
        diagonal - boolean indicating whether to plot y = x line
        hline - value to add horziontal line at
        vline - value to add vertical line at
        title - title to add to scatter
        legend_out - boolean indicating whether legend should be placed outside axis
        hide_legend - boolean indicating if legend should be hidden
        text_kws - keyword args for text
        arrowprops - properties to apply to arrows for labels
        return_ax - boolean indicating whether to return the axis
        ax - matplotlib axis to add the scatter to
        use_lims - x and y lims to apply
        adj_kws - adjust text keyword args
        kwargs - other scatterplot keyword args
    '''
    #format labels
    if isinstance(labels, list) or isinstance(labels, pd.Index):
        labels = pd.Series(labels, labels)
    #align x,y, and hue into dataframe
    data = pd.concat([x, y, hue, labels.rename("Label") if labels is not None else labels], axis=1)
    assert data.columns.is_unique, "non-unique names for inputs"
    #setup new figure if ax not provided
    if ax is None:
        plt.figure(figsize=figsize)
        plt.title(title)
        ax = plt.gca()
    #align hue_order and style_order if necessary
    if "hue_order" in kwargs and style_hue:
        kwargs["style_order"] = kwargs["hue_order"]
    #create scatter
    sns.scatterplot(
        data=data, x=x.name, y=y.name, 
        hue=hue.name if hue is not None else None, 
        style=hue.name if style_hue else None, 
        ax=ax, 
        **kwargs
    )
    #add a diagonal y=x line
    if diagonal:
        clean_data = data[[x.name, y.name]].dropna()
        lims = (clean_data.min().min(),clean_data.max().max())
        ax.plot(lims, lims, color=ACCENT, alpha=0.5, linestyle="--", zorder=-1)
    #add horizontal/vertical lines
    if hline is not None:
        axline(ax, "h", hline)
    if vline is not None:
        axline(ax, "v", vline)
    #handle legend appropriately
    handles, labels = plt.gca().get_legend_handles_labels()
    if hue is not None:
        if not hide_legend:
            if legend_out:
                legend_out = dict(bbox_to_anchor=(1.04,1), loc="upper left")
            else:
                legend_out = dict()
        plt.gca().legend(title=hue.name, **legend_out)
    #set x/y li,s
    if use_lims is not None:
        plt.xlim(use_lims[0])
        plt.ylim(use_lims[1])
    #add labels
    if "Label" in data.columns:
        from adjustText import adjust_text
        texts = [
            ax.text(data.loc[i, x.name], data.loc[i, y.name], data.loc[i, "Label"], **text_kws)
            for i in data.index[data["Label"].notnull()]
        ]
        adjust_text(texts, arrowprops=arrowprops, ax=ax, **adj_kws)
    #optionally return the axis
    if return_ax:
        return plt.gca()

def sem(s):
    if len(s)<3:
        return [None, None]
    else:
        avg = s.mean()
        sem = scipy.stats.sem(s)
        return [avg-sem, avg+sem]

def fancybar(
    data, x, y, hue=None, palette=None, hue_order=None, order=None, ax=None, capsize=0.1, edgewidth=0.1, s=5, 
    errorbar="se", dodge=True, width=0.75, err_kws=dict(linewidth=0.5), bar_alpha=0.5
):
    '''
    Produce a nicely styled barplot with points.

    Parameters:
        data - pandas DataFrame to plot
        x - name of column in data to use for x-axis
        y - name of column in data to use for y-axis
        hue - name of column in data to color by
        hue_order - order of colors
        order - order of values on axis
        ax - matplotlib axis to add the scatter to
        capsize - size of errorbar caps
        edgewidth - size of errorbar lines
        s - size of points
        errorbar - type of errorbar
        dodge - whether to dodge colors as separate columns
        width - width of bars
        err_kws - other keywords for errorbars
        bar_alpha - transperancy of errorbar
    '''
    #create a stripplot for points
    curr_ax = sns.stripplot(
        data=data, x=x, y=y, hue=hue, dodge=0.1 if dodge else False, jitter=True, order=order,
        hue_order=hue_order, palette=palette, s=s, alpha=0.75, linewidth=0.5, edgecolor=ACCENT,
        zorder=1, ax=ax
    )
    #create bars if hue order is provided, to hide the colors from the legend
    if hue_order is not None:
        sns.barplot(
            data=data.replace({x:"_"+x for x in hue_order}), order=order,
            x=x, y=y, hue=hue, dodge=0.5 if dodge else False, hue_order=["_"+x for x in hue_order], palette=palette,
            alpha=bar_alpha, width=width, capsize=capsize, err_kws=err_kws,
            errorbar=errorbar, linestyle="none", zorder=0, ax=ax
        )
    #otherwise just add bar
    else:
        sns.barplot(
            data=data, order=order,
            x=x, y=y, hue=hue, dodge=0.5 if dodge else False, palette=palette,
            alpha=bar_alpha, width=width, capsize=capsize, err_kws=err_kws,
            errorbar=errorbar, linestyle="none", zorder=0, ax=ax
        )
    return curr_ax

def fan_plot(
    data, y, groups=None, hue=None, colors=None, palette="tab20",
    inverse=False, log=True, gap=1, label=True, label_dist=1.5, 
    h_lines={0:0}, h_line_angle=0.1, 
    figsize=(5, 5), title=None, fontsize=6, ax=None, labelpad=30
):
    '''
    Plots a circular barplot, aka a fan plot.  

    Parameters:
        data - pandas DataFrame containing y and hue columns
        y - column to use for y-axis (heights)
        groups - column to group on
        hue - column to color by
        colors - colors to assign to categories
        palette - color palette to use
        inverse - whether to inverse the heights by negating them
        log - whether to rescale axis
        gap - how large of a gap to leave in center of plot
        label - whether to label groups
        label_dist - how far from axis to put labels
        h_lines - horizontal lines to plot (dict with labels as key and height as value)
        h_line_angle - angle to display horizontal line labels on
        figsize - size for figure
        title - title for figure
        fontsize - size of text for labels
        ax - axis to plot on
        labelpad - how much padding to add to labels
    '''
    #order data
    data = data.sort_values([groups, y]) if groups is not None else data.sort_values(y)

    #values for the x and y axis
    angles = np.linspace(0, 2 * np.pi, len(data), endpoint=False)
    heights = data[y].values * (-1 if inverse else 1)
    
    #initialize layout in polar coordinates
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    #pick colors
    if hue is not None:
        if colors is None:
            categories = pd.Categorical(data[hue]).codes
            colormap = sns.color_palette(palette, len(categories))
            colors = np.array(colormap)[categories]
        else:
            colors = data.merge(
                colors.to_frame("Colors"), how="left", left_on=hue, right_index=True
            )["Colors"].values
    
    #setup polar axis
    if log:
        ax.set_rscale('symlog')
    ax.set_theta_offset(np.pi / 2)
    
    #add lines
    ax.vlines(angles, 0 + gap, heights + gap,  lw=0.9, color=colors)
    
    #remove extra plot features
    ax.spines["start"].set_color("none")
    ax.spines["polar"].set_color("none")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticklabels([])

    #add horizontal lines
    h_line_angles = np.linspace(0, 2 * np.pi, 200)
    for label, line in h_lines.items():
        ax.plot(h_line_angles, np.repeat(line + gap, 200), color="black", lw=0.7)
        ax.text(
            h_line_angle, line+(gap*0.80), label, color="black", fontsize=fontsize+2
        )

    #label groups
    if groups is not None:
        group_size = data.groupby(groups).size()
        start = 0 
        for group, size in group_size.items():
            #get middle of group
            middle = int((size / 2) + start)
            #get rotation and alignment
            rotation = np.rad2deg(angles[middle]+np.pi / 2)
            if angles[middle] <= np.pi:
                alignment = "right"
                rotation = rotation + 180
            else: 
                alignment = "left"
            #add label
            if label:
                ax.text(
                    angles[middle], gap*label_dist, group, fontsize=fontsize,
                    ha=alignment, va="center", rotation=rotation,
                    rotation_mode="anchor"
                )
            start += size

    ax.set_xlabel(title, labelpad=labelpad)


## Western blots ##

def crop_and_box_wb(img_path, crops, cmap=None, save=False):
    '''
    Crop a western blot image, adding box to show which part has been cropped
    '''
    #if image is found, proceed
    if os.path.exists(img_path):
        #get image and plot it
        img = np.asarray(Image.open(img_path).convert('RGB'))
        fig = plt.figure()
        plt.gca().imshow(img, aspect='equal', cmap=cmap)
        cropped = {}
        #show image without cropping
        if crops is None:
            plt.show()
            return None
        #get crops and add corresponding boxes
        for label, c in crops.items():
            cropped[label] = img[c[0]:c[1], c[2]:c[3]]
            plt.gca().add_patch(plt.Rectangle((c[2], c[0]), c[3]-c[2], c[1]-c[0], ls="--", ec="r", fc="none"))
            plt.gca().annotate(label, (c[2], c[0]), (0, 10), textcoords='offset pixels', color="r")
        #save boxed figure
        if save:
            new_path = img_path.split(".")
            plt.savefig(new_path[0]+"_boxed")
        plt.gca().set_axis_off()
        plt.show()
    #otherwise create empty crops of same size
    elif crops is not None:
        cropped = {label:np.ones((c[1]-c[0], c[3]-c[2]))*255 for label, c in crops.items()}
    else:
        cropped = None
    return cropped


LABEL_AX_SIZE = (1, 5)
def get_wb_height_ratios(imgs, add_label_ratio=LABEL_AX_SIZE):
    '''
    Determine heights of western blot panels to keep images proportional.
    '''
    shapes = [add_label_ratio]+[img.shape for img in imgs]
    max_width = pd.Series([shape[1] for shape in shapes]).max()
    max_height = pd.Series([shape[0] for shape in shapes]).max()
    return [(shape[0] / max_height) / (shape[1] / max_width) for shape in shapes]

def add_wb_labels(
    ax, levels, line_gap=0.02, level_params={"DOX":{"fontsize":7, "ha":"center", "va":"center"}}, 
    linewidth=0.55, size=LABEL_AX_SIZE, groups=None
):
    '''
    Create an axis with WB annotations.
    '''
    #set aspect to make width match with ratio
    ax.set_axis_off()
    ax.set_aspect("equal")
    ax.set_ylim(0, size[0])
    ax.set_xlim(0, size[1])
    #determine size of gap between lines
    line_gap = size[1]*line_gap
    #add group labels
    if groups is not None:
        #set distance between axes 
        dy = size[0] / (len(levels)+2)
        i=2
        for (group, (x0,x1)) in groups:
            #add text and underline for each group
            params = {"fontsize":5, "ha":"center", "va":"center"} if "groups" not in level_params.keys() else level_params["groups"]
            left, right = x0*size[1], x1*size[1]
            middle = left+((right-left)/2)
            ax.annotate(
                group, (middle, 1),  weight="bold", annotation_clip=False, 
                **params
            )
            sns.lineplot(
                x=[left+line_gap, right-line_gap], y=[1-dy, 1-dy], ax=ax, color="k", 
                linewidth=linewidth*1.5
            )
    else:
        #set distance between axes 
        dy = size[0] / len(levels)
        i = 0
    #for each level of labels
    for level, labels in levels:
        #get parameters and height
        params = {"fontsize":5, "ha":"center", "va":"center"} if level not in level_params.keys() else level_params[level]
        height = 1-(i*dy)
        #add the y-axis label
        if level != "_":
            ax.annotate(
                level, (0, height), annotation_clip=False, fontsize=5, ha="right", va="center"
            )
        #figure out spacing of annots/lines, then add them
        n = size[1]/len(labels)
        for j, label in enumerate(labels):
            left, middle, right = (j*n), ((j+0.5)*n), ((j+1)*n)
            if level != "_":
                ax.annotate(
                    label, (middle, height), annotation_clip=False, **params
                )
            else:
                sns.lineplot(
                    x=[left+line_gap, right-line_gap], y=[height, height], ax=ax, color="k", 
                    linewidth=linewidth
                )
        i+=1

def wb_image(ax, label, img, ladder=None, fontsize=5):
    '''
    Prepare a WB image with ladder if desired.
    '''
    #label image, then display it with no ticks/tick labels
    ax.set_ylabel(label, fontsize=fontsize, rotation=0, va="center", ha="right")
    ax.imshow(img, aspect="equal")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    for side in 'left bottom top right'.split():
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.75)
    #add ladder
    if ladder is not None:
        ax2 = ax.secondary_yaxis('right')
        ax2.tick_params('both', length=3, width=1, which='major', pad=1)
        ax2.set_yticks([(1-p)*img.shape[0] for p in ladder.values()])
        ax2.set_yticklabels(list(ladder.keys()), fontsize=fontsize)
        ax2.spines["right"].set_visible(False)

