import plotly.express as px

from .utils import _default_to_if_None_and_logic
from .colors import colors
from ..utils import do_call
from ..list import unique, order

def scatter_plotly(
    data_frame, x_by: str, y_by: str, color_by: str,
    px_args: dict = {},
    size = 5, color_panel: list = colors,
    color_order_legend: str = 'increasing',
    plot_title: str = None, legend_title: str = None,
    xlab: str = None, ylab: str = None
    ):
    """
    Produces a scatter plot using plotly.express.scatter based on the pandas 'data_frame' given.
    'x_by', 'y_by', and 'color_by' should indicate columns of 'data_frame' to use for x/y/color.
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.scatter' call.
    'size' sets the size of points.  Can be either a number directly or the name of a column of 'data_frame'.
    'color_panel' (string list) sets the colors when 'color_by' references discrete data
    'color_order_legend' ('increasing', 'decreasing', or 'unordered') sets the ordering of keys in the legend, when 'color_by' references discrete data
    'plot_title', 'legend_title', 'xlab', and 'ylab' set titles.
    """

    # Parse dependent defaults (given 'None' b/c python)
    xlab = _default_to_if_None_and_logic(xlab, x_by)
    ylab = _default_to_if_None_and_logic(ylab, y_by)
    plot_title = _default_to_if_None_and_logic(plot_title, color_by)
    legend_title = _default_to_if_None_and_logic(legend_title, color_by)

    # Add to px_args (Can probably remove this section once I learn python syntax better!)
    px_args['data_frame'] = data_frame
    px_args['x'] = x_by
    px_args['y'] = y_by
    px_args['color'] = color_by
    px_args['color_discrete_sequence'] = color_panel

    # Set legend key order.
    # FUTURE! Make this work for bool too.  And maybe other types
    if not all(map(lambda x: isinstance(x, str), data_frame[color_by])):
        color_order_legend = 'unordered'
    if (color_order_legend == 'increasing' or 'decreasing'):
        categories = order(unique(data_frame[color_by]))
        px_args['category_orders'] = { color_by: categories }
    if (color_order_legend == 'decreasing'):
        px_args['category_orders'] = { color_by: list(reversed(categories)) }
    
    # Make plot
    fig = do_call(px.scatter, px_args)

    fig.update_layout(
        title_text=plot_title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        legend_title=legend_title,
        legend= {'itemsizing': 'constant'}
    )

    # Tweaks
    # fig.update_coloraxes(colorbar_title_text=color_by)
    fig.update_traces(marker={'size': size})

    return fig

