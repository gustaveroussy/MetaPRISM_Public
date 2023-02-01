import pandas as pd

# plotly
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def get_table_counts(df, x, add_pct=False):
    """Prepare table of counts of values in column 'x'.

    Parameters
    ----------
    df : dataframe
    x: str
        Name of the tumor type column.
    add_pct: bool, optional

    Returns
    -------
    df_count: dataframe
    """
    df_count = df[x].value_counts()
    df_count = df_count.to_frame("Count").reset_index()
    df_count.columns = [x, "Count"]
    if add_pct:
        df_count["Percentage"] = df_count["Count"]/df_count["Count"].sum()
    return df_count



def select_tumor_types(df_cln, min_count=None, tumor_types_keep=None, tumor_types_drop=None, x="Project_TCGA_More"):
    """Select tumor types from clinical table.

    Parameters
    ----------
    df_cln : dataframe
    min_count : int, optional
    tumor_types_drop : list-like, optional
    tumor_types_keep : list-like, optional
    x: str, optional
        Name of the tumor type column.

    Returns
    -------
    df_cln_sub: dataframe
    """
    if tumor_types_drop is not None:
        df_cln = df_cln.loc[~df_cln[x].isin(tumor_types_drop)]

    if tumor_types_keep is None:
        df_counts = get_table_counts(df_cln, x)
        tumor_types_keep = df_counts.loc[df_counts["Count"]>=min_count, x].tolist()

    df_cln_sub = df_cln.loc[df_cln[x].isin(tumor_types_keep)]
    return df_cln_sub


# DONUT PLOT ===========================================================================================================

def plot_donut_count(df_counts, name="", field_label="Label", field_count="Count", field_color="Color", hole=0.4,
                     size=12, **kwargs):
    df_counts = df_counts.sort_values(by=[field_label])
    pull = df_counts["Pull"] if "Pull" in df_counts.columns else None
    pie = go.Pie(labels=df_counts[field_label],
                 values=df_counts[field_count],
                 insidetextfont=dict(family="Helvetica", size=size),
                 outsidetextfont=dict(family="Helvetica", color="black", size=size),
                 pull=pull,
                 textinfo="label",
                 marker=dict(line=dict(color='#000000', width=0.35), colors=df_counts[field_color]),
                 sort=False, hole=hole, name="%s<br>%s" % (name.upper(), f'{df_counts[field_count].sum():,}'),
                 **kwargs)
    return pie


def annotation_centered_name(fig_data, size=20):
    x_mid = (fig_data["domain"]["x"][0] + fig_data["domain"]["x"][1])/2
    y_mid = (fig_data["domain"]["y"][0] + fig_data["domain"]["y"][1])/2
    annotation = dict(text=fig_data["name"], x=x_mid, y=y_mid, showarrow=False,
                      font=dict(size=size, family="Helvetica", color="black"), xanchor="center", yanchor="middle")
    return annotation


def draw_plot_donut_count(cohort, df_counts, size=12, size_annotation=20, **kwargs):
    fig = make_subplots(1, 1, specs=[[{'type':'domain'}]])
    trace = plot_donut_count(df_counts, name=cohort, size=size, **kwargs)
    fig.add_trace(trace, 1, 1)

    annotation = annotation_centered_name(fig_data=fig.to_dict()["data"][0], size=size_annotation)
    margin = dict(t=0, l=0, b=0, r=0)
    fig.update_layout(annotations=[annotation], showlegend=False, margin=margin)

    return fig
