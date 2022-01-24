import plotly.express as px

from .utils import _default_to_if_make_and_logic
from .colors import colors
from ..list import unique, order

def scatter_plotly(
    data_frame, x_by: str, y_by: str,
    color_by: str = "make",
    px_args: dict = {},
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
    """

    # Parse dependent defaults
    xlab = _default_to_if_make_and_logic(xlab, x_by)
    ylab = _default_to_if_make_and_logic(ylab, y_by)
    plot_title = _default_to_if_make_and_logic(plot_title, color_by)
    legend_title = _default_to_if_make_and_logic(legend_title, color_by)

    # Add to px_args
    px_args['data_frame'] = data_frame
    px_args['x'] = x_by
    px_args['y'] = y_by
    px_args['color_discrete_sequence'] = color_panel
    px_args['hover_data'] = hover_data
    
    # Set coloring if given data to color by
    if color_by!="make":
        px_args['color'] = color_by
        # Also set plotting/legend key order
        discrete_color = any(map(lambda x: isinstance(x, (str, bool)), data_frame[color_by]))
        if (color_order == 'increasing' or color_order == 'decreasing'):
            categories = order(unique(data_frame[color_by]))
            if (color_order == 'increasing'):
                px_args['category_orders'] = { color_by: categories }
                if order_when_continuous_color and not discrete_color:
                    px_args['data_frame'] = data_frame.iloc[ order(data_frame[color_by], return_indexes=True) ]
            else:
                px_args['category_orders'] = { color_by: list(reversed( categories )) }
                if order_when_continuous_color and not discrete_color:
                    px_args['data_frame'] = data_frame.iloc[ list(reversed( order(data_frame[color_by], return_indexes=True) )) ]
    else:
        plot_title = _default_to_if_make_and_logic(plot_title, "")
    
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

