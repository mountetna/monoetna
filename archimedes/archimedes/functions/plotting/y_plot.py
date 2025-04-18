from typing import List, Union
import plotly.express as px
import pandas as pd

from .utils import _leave_default_or_none, _which_rows
from .colors import colors
from ..list import order, unique

def y_plotly(
    data_frame: pd.DataFrame,
    x_by: str,
    y_by: str,
    color_by: str = "make",
    plots: List[str] = ["violin", "box"],
    px_args: dict = {},
    rows_use = None,
    x_order: Union[str, list] = 'unordered',
    y_scale = "linear",
    color_panel: list = colors,
    xlab: str = "make",
    ylab: str = "make",
    plot_title: str = "make",
    legend_title: str = "make"):
    """
    Produces violin and/or box plots using plotly.express.violin, or plotly.express.box, based on the pandas 'data_frame' given.
    'x_by', 'y_by' should indicate columns of 'data_frame' to use for x/y axes data. 'x_by' should denote discrete groupings, while 'y-by' should denote numerical data.
    'color_by' can indicate a column of 'data_frame' to use for establishing fill colors of violin/boxplots and can be used to create subgroupings, or highlight supergroupings, of 'x_by' groups. By default, if left as "make", it will be passed the 'x_by' data.
    'plots' sets which data representations to use, and should be a list containing "violin", "box", or both of these.
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.violin' or 'plotly.express.box' call.
    'color_panel' (string list) sets the colors to use for violin/boxplot fills.
    'plot_title', 'legend_title', 'xlab', and 'ylab' set titles.
    'rows_use'
    'x_order', either "unordered", "increasing", "decreasing", or a list[str] with the desired order for x-axis groups. 
    'y_scale', String, 'linear', 'log10', or 'log10(val+1)'. Controls whether this axes should be log scaled, and if so, whether 1 should be added to all values first in order to let zeros be okay to plot.
    """
    
    # Parse dependent defaults
    color_by = _leave_default_or_none(color_by, x_by)
    xlab = _leave_default_or_none(xlab, x_by)
    ylab = _leave_default_or_none(ylab, y_by)
    plot_title = _leave_default_or_none(plot_title, ylab)
    legend_title = _leave_default_or_none(legend_title, color_by)
    
    # data_frame edits
    df = data_frame.copy()
    rows_use = _which_rows(rows_use, data_frame)
    df = df.loc[rows_use]
    if y_scale=="log10(val+1)":
        y_scale="log10"
        df[y_by] += 1
    
    # Add to px_args and convert from our variables names to px.violin/bar variables names. 
    px_args["data_frame"] = df
    px_args["x"] = x_by
    px_args["y"] = y_by
    px_args["color"] = color_by
    px_args["color_discrete_sequence"] = color_panel
    px_args["log_y"] = y_scale=="log10"
    if x_order!="unordered":
        categories = {}
        if isinstance(x_order, list):
            categories[x_by] = x_order
        elif x_order in ["increasing","decreasing"]:
            categories[x_by] = order(unique(df[x_by]), decreasing=x_order=="decreasing")
        px_args["category_orders"]=categories
    
    # Make Plot
    if "violin" in plots:
        #use conditional here for checking if box is in plots
        px_args["box"]= "box" in plots
        fig = px.violin(
            **px_args
            )
    elif "box" in plots:
        fig = px.box(
            **px_args
            )
    
    # Tweaks
    fig.update_layout(
        title_text=plot_title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        legend_title_text=legend_title,
        legend= {'itemsizing': 'constant'}
    )
    return fig

