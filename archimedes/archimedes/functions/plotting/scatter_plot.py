import plotly.express as px

from plotnine import ggplot, aes, theme_bw, ggtitle, xlab, ylab, facet_wrap, facet_grid, scale_y_continuous, scale_x_continuous, scale_colour_gradient, scale_colour_manual, guides, guide_legend, geom_point, geom_density_2d, scale_shape_manual

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
    x_title: str = "make", y_title: str = "make",
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
    'plot_title', 'legend_title', 'x_title', and 'y_title' set titles.
    """

    # Parse dependent defaults
    x_title = _default_to_if_make_and_logic(x_title, x_by)
    y_title = _default_to_if_make_and_logic(y_title, y_by)
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
        xaxis_title=x_title,
        yaxis_title=y_title,
        legend= {'itemsizing': 'constant'}
    )

    # Tweaks
    fig.update_coloraxes(colorbar_title_text=legend_title)
    fig.update_traces(marker={'size': size}, )

    return fig

def _add_contours(
    fig, data, x_by, y_by, color, linetype = 1):
    # Add contours based on the density of cells/samples
    # (Dim and Scatter plots)
    
    return fig + geom_density_2d(
        data = data,
        mapping = aes(x = x_by, y = y_by),
        color = color,
        linetype = linetype,
        na_rm = True)

def scatter_plotnine(
    data_frame, x_by: str, y_by: str,
    color_by: str = '',
    size = 5, color_panel: list = colors,
    color_order: str = 'increasing',
    plot_order: str = 'increasing',
    plot_title: str = "make", legend_title: str = "make",
    x_title: str = "make", y_title: str = "make",
    plot_theme = theme_bw(),
    opacity = 1,
    min_color = "#F0E442",
    max_color = "#0072B2",
    min_value = None,
    max_value = None,
    legend_color_breaks = "make",
    legend_color_breaks_labels = "make",
    legend_color_size = 5,
    split_by = [],
    do_contour = False,
    contour_color = "black",
    contour_linetype = 'solid',
    shape_by = None,
    shape_panel = ['o','s','^','D','v','*'], 
    legend_shape_title = "make",
    legend_shape_size = 5,
    y_scale = scale_y_continuous,
    x_scale = scale_x_continuous
    ):
    
    # Parse dependent defaults
    x_title = _default_to_if_make_and_logic(x_title, x_by)
    y_title = _default_to_if_make_and_logic(y_title, y_by)
    plot_title = _default_to_if_make_and_logic(plot_title, color_by)
    legend_title = _default_to_if_make_and_logic(legend_title, color_by)
    legend_shape_title = _default_to_if_make_and_logic(legend_shape_title, shape_by)
    
    ### Start plot, with data and theming
    fig = (ggplot(data_frame) +
        ylab(y_title) +
        xlab(x_title) +
        plot_theme +
        x_scale() +
        y_scale()
    )
    
    if plot_title!=None:
        fig += ggtitle(plot_title)
    
    aes_args = {'x': x_by, 'y': y_by}
    geom_args = {
        'data': data_frame,
        'size': size,
        'alpha': opacity}
    
    if color_by!='':
        aes_args['color'] = color_by
        
        if isinstance(data_frame[color_by][0], (int, float, complex)):
            scale_args= {
                'name': legend_title,
                'low': min_color,
                'high': max_color,
                'limits': (
                    [min_value,min(data_frame[color_by])][min_value==None]
                    [max_value,max(data_frame[color_by])][max_value==None])
            }
            if legend_color_breaks!="make":
                scale_args['breaks'] = legend_color_breaks,
            if legend_color_breaks_labels!="make":
                scale_args['labels'] = legend_color_breaks_labels
            fig += scale_colour_gradient(**scale_args)
        else:
            fig += scale_colour_manual(
                name = legend_title,
                values = color_panel)
            fig += guides(color = guide_legend(override_aes = {'size':legend_color_size}))
    
    if shape_by!=None:
        aes_args['shape'] = shape_by
        fig += scale_shape_manual(
                values = shape_panel,
                name = legend_shape_title)
        fig += guides(shape = guide_legend(override_aes = {'size':legend_shape_size}))
    else:
        geom_args['shape'] = shape_panel[0]
    
    ### Add Data
    geom_args['mapping'] = aes(**aes_args)
    fig += geom_point(**geom_args)
    
    ### Extra tweaks
    # Faceting
    if len(split_by)==1:
        fig += facet_wrap(split_by)
    if len(split_by)==2:
        fig += facet_grid(split_by)
    # Contours
    if do_contour:
        fig = _add_contours(fig, data_frame, x_by, y_by, contour_color, contour_linetype)
    
    return fig
    
    