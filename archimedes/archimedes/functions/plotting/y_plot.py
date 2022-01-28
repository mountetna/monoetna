import plotly.express as px

from .utils import _leave_default_or_null, _which_rows
from .colors import colors

def y_plotly(
    data_frame,
    x_by,
    y_by,
    color_by = "make",
    plots = ["violin", "box"],
    px_args: dict = {},
    rows_use = None,
    y_scale = "as is",
    color_panel: list = colors,
    xlab="make",
    ylab="make",
    plot_title="make",
    legend_title ="make"):
    """
    Produces violin and/or box plots using plotly.express.violin, or plotly.express.box, based on the pandas 'data_frame' given.
    'x_by', 'y_by' should indicate columns of 'data_frame' to use for x/y axes data. 'x_by' should denote discrete groupings, while 'y-by' should denote numerical data.
    'color_by' can indicate a column of 'data_frame' to use for establishing fill colors of violin/boxplots and can be used to create subgroupings, or highlight supergroupings, of 'x_by' groups. By default, if left as "make", it will be passed the 'x_by' data.
    'plots' sets which data representations to use, and should be a list containing "violin", "box", or both of these.
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.violin' or 'plotly.express.box' call.
    'color_panel' (string list) sets the colors to use for violin/boxplot fills.
    'plot_title', 'legend_title', 'xlab', and 'ylab' set titles.
    'rows_use'
    'y_scale', String, 'as is' or 'log10'. Controls whether this axes should be log scaled. (Not coded as boolean in anticipation of a static plotter system offering extended options)
    """
    
    # Parse dependent defaults
    color_by = _leave_default_or_null(color_by, x_by)
    xlab = _leave_default_or_null(xlab, x_by)
    ylab = _leave_default_or_null(ylab, y_by)
    plot_title = _leave_default_or_null(plot_title, ylab)
    legend_title = _leave_default_or_null(legend_title, color_by)
    
    # data_frame edits
    df = data_frame.copy()
    rows_use = _which_rows(rows_use, data_frame)
    df = df.loc[rows_use]
    
    # Add to px_args and convert from our variables names to px.violin/bar variables names. 
    px_args["data_frame"] = df
    px_args["x"] = x_by
    px_args["y"] = y_by
    px_args["color"] = color_by
    px_args["color_discrete_sequence"] = color_panel
    px_args["log_y"] = y_scale=="log10"
    
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

