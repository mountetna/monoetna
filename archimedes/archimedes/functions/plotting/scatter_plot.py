import plotly.express as px

from .utils import _leave_default_or_none, _which_rows, _is_discrete
from .colors import colors
from ..list import unique, order

def scatter_plotly(
    data_frame, x_by: str, y_by: str,
    color_by: str = "make",
    px_args: dict = {},
    rows_use = None,
    x_scale = "as is", y_scale = "as is",
    size = 5, color_panel: list = colors,
    color_order: str = 'increasing',
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
    'color_order' ('increasing', 'decreasing', or 'unordered') sets the ordering of keys in the legend, when 'color_by' references discrete data
    'plot_title', 'legend_title', 'xlab', and 'ylab' set titles.
    'rows_use',
    'x_scale', 'y_scale'. String, 'as is', 'log10', or 'log10(val+1)'. Controls whether these axes should be log scaled, and if so, whether 1 should be added to all values first in order to let zeros be okay to plot.
    """

    # Parse dependent defaults
    xlab = _leave_default_or_none(xlab, x_by)
    ylab = _leave_default_or_none(ylab, y_by)
    plot_title = _leave_default_or_none(plot_title, color_by)
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
        discrete_color = _is_discrete(df[color_by])
        if (color_order == 'increasing' or color_order == 'decreasing'):
            categories = order(unique(df[color_by]))
            if (color_order == 'increasing'):
                px_args['category_orders'] = { color_by: categories }
                if order_when_continuous_color and not discrete_color:
                    px_args['data_frame'] = df.iloc[ order(df[color_by], return_indexes=True) ]
            else:
                px_args['category_orders'] = { color_by: list(reversed( categories )) }
                if order_when_continuous_color and not discrete_color:
                    px_args['data_frame'] = df.iloc[ list(reversed( order(df[color_by], return_indexes=True) )) ]
    else:
        plot_title = _leave_default_or_none(plot_title, "")
    
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

