from typing import Union
import pandas as pd
import plotly.express as px

from .utils import _leave_default_or_none, _which_rows, _is_discrete
from .colors import colors
from ..list import unique, order

def scatter_plotly(
    data_frame: pd.DataFrame, x_by: str, y_by: str,
    color_by: str = "make",
    px_args: dict = {},
    rows_use = None,
    x_scale = "linear", y_scale = "linear",
    size = 5, color_panel: list = colors,
    color_order: Union[str, list] = 'increasing',
    order_when_continuous_color: bool = False,
    plot_title: str = "make", legend_title: str = "make",
    xlab: str = "make", ylab: str = "make",
    hover_data: str = None
    ):
    """
    Produces a scatter plot using plotly.express.scatter based on the pandas 'data_frame' given.
    'x_by', 'y_by' should indicate columns of 'data_frame' to use for x/y axes data.
    'color_by' should indicate a column of 'data_frame' to use for coloring the data points, OR, default, if left as "make" data points will all be a single color.
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.scatter' call.
    'size' sets the size of points.  Can be either a number directly or the name of a column of 'data_frame'.
    'color_panel' (string list) sets the colors when 'color_by' references discrete data
    'order_when_continuous_color'  sets the ordering of data points from back to front.
    'color_order' ('increasing', 'decreasing', 'unordered', or a list) sets the render ordering from back to front as well as, when 'color_by' references discrete data, the ordering of colors & keys in the legend. When a list, should be the data values in their desired order.
    'plot_title', 'legend_title', 'xlab', and 'ylab' set titles.
    'rows_use',
    'x_scale', 'y_scale'. String, 'linear', 'log10', or 'log10(val+1)'. Controls whether these axes should be log scaled, and if so, whether 1 should be added to all values first in order to let zeros be okay to plot.
    """

    # Parse dependent defaults
    xlab = _leave_default_or_none(xlab, x_by)
    ylab = _leave_default_or_none(ylab, y_by)
    plot_title = _leave_default_or_none(plot_title, "" if color_by=="make" else color_by)
    legend_title = _leave_default_or_none(legend_title, color_by)
    
    # data_frame edits
    df = data_frame.copy()
    rows_use = _which_rows(rows_use, df)
    df = df.loc[rows_use]
    if x_scale=="log10(val+1)":
        df[x_by] += 1
        x_scale="log10"
    if y_scale=="log10(val+1)":
        df[y_by] += 1
        y_scale="log10"

    # Add to px_args
    px_args['data_frame'] = df
    px_args['x'] = x_by
    px_args['y'] = y_by
    px_args['color_discrete_sequence'] = color_panel
    px_args['hover_data'] = hover_data
    px_args["log_x"] = x_scale=="log10"
    px_args["log_y"] = y_scale=="log10"
    
    # Set coloring if given data to color by
    if color_by!="make":
        px_args['color'] = color_by
        # Also set plotting/legend key order
        categories = {}
        if isinstance(color_order, list):
            categories[color_by] = color_order
        elif (color_order in ['increasing','decreasing']):
            categories[color_by] = order(unique(df[color_by]), decreasing=color_order=="decreasing")
            if order_when_continuous_color and not _is_discrete(df[color_by]):
                px_args['data_frame'] = df.iloc[ order(df[color_by], return_indexes=True, decreasing=color_order=="decreasing") ]
        px_args['category_orders'] = categories
    
    # Make plot
    fig = px.scatter(**px_args)

    fig.update_layout(
        title_text=plot_title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        legend= {'itemsizing': 'constant'}
    )

    # Tweaks
    fig.update_coloraxes(colorbar_title_text=legend_title)
    fig.update_traces(marker={'size': size}, )

    return fig

