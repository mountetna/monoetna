import plotly.express as px

from .utils import _default_to_if_make_and_logic
from .colors import colors

def y_plotly(
    data_frame,
    x_by,
    y_by,
    color_by = "make",
    plots = ["violin", "box"],
    px_args: dict = {},
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
    """
    
    # Parse dependent defaults
    color_by = _default_to_if_make_and_logic(color_by, x_by)
    xlab = _default_to_if_make_and_logic(xlab, x_by)
    ylab = _default_to_if_make_and_logic(ylab, y_by)
    plot_title = _default_to_if_make_and_logic(plot_title, ylab)
    legend_title = _default_to_if_make_and_logic(legend_title, color_by)
    
    # Add to px_args and convert from our variables names to px.violin/bar variables names. 
    px_args["data_frame"] = data_frame
    px_args["x"] = x_by
    px_args["y"] = y_by
    px_args["color"] = color_by
    px_args["color_discrete_sequence"] = color_panel
    
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

