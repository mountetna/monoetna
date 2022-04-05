import plotly.express as px
import numpy as np

from plotnine import ggplot, aes, theme_bw, ggtitle, xlab, ylab, scale_y_continuous, scale_x_continuous, scale_colour_gradient, scale_colour_manual, guides, guide_legend, geom_point, scale_shape_manual, theme, element_blank

from .utils import _leave_default_or_none, _is_discrete, _is_continuous, _is_logical, _is_integer, _scale, _all_rows, _which_rows
from .utils_col_getters import _isCol, _col, _colLevels
from .utils_plot_mods import _add_contours, _add_splitting, _remove_legend
from .colors import colors
from ..list import unique, order

def scatter_plotly(
    data_frame, x_by: str, y_by: str,
    color_by: str = "make",
    px_args: dict = {},
    rows_use = None,
    x_scale = "linear", y_scale = "linear",
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
    x_title = _leave_default_or_none(x_title, x_by)
    y_title = _leave_default_or_none(y_title, y_by)
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
        xaxis_title=x_title,
        yaxis_title=y_title,
        legend= {'itemsizing': 'constant'}
    )

    # Tweaks
    fig.update_coloraxes(colorbar_title_text=legend_title)
    fig.update_traces(marker={'size': size}, )

    return fig

def scatter_plotnine(
    data_frame,
    x_by: str,
    y_by: str,
    color_by: str = '',
    shape_by: str = '',
    split_by = [],
    size = 1,
    rows_use = None,
    show_others = True,
    color_adjustment = None,
    color_adj_fxn = None,
    split_show_all_others = True,
    opacity: float = 1,
    color_panel: list = colors,
    colors = None,
    split_nrow = None,
    split_ncol = None,
    split_adjust = {},
    shape_panel = ['o','s','^','D','v','*'],
    # rename_color_groups = None,
    # rename_shape_groups = None,
    min_color = "#F0E442",
    max_color = "#0072B2",
    min_value = None,
    max_value = None,
    plot_order: str = 'unordered',
    x_title: str = "make",
    y_title: str = "make",
    plot_title: str = "make",
    plot_theme = theme_bw(),
    do_contour = False,
    contour_color = "black",
    contour_linetype = 'solid',
    # add_trajectory_by_groups = None,
    # add_trajectory_curves = None,
    # trajectory_group_by = None,
    # trajectory_arrow_size = 0.15,
    # do_letter = False,
    # do_ellipse = False,
    # do_label = False,
    # labels_size = 5,
    # labels_highlight = True,
    # labels_repel = True,
    # labels_split_by = "make",
    legend_show = True,
    legend_color_title: str = "make",
    legend_color_size = 5,
    legend_color_breaks = "make",
    legend_color_breaks_labels = "make",
    legend_shape_title = "make",
    legend_shape_size = 5,
    show_grid_lines = True,
    y_scale = scale_y_continuous,
    x_scale = scale_x_continuous,
    data_out = False
    ):
    
    df = data_frame.copy()
    
    # Parse dependent defaults
    x_title = _leave_default_or_none(x_title, x_by)
    y_title = _leave_default_or_none(y_title, y_by)
    plot_title = _leave_default_or_none(plot_title, color_by)
    legend_color_title = _leave_default_or_none(legend_color_title, color_by)
    legend_shape_title = _leave_default_or_none(legend_shape_title, shape_by)
    if colors!= None:
        color_panel = color_panel[colors]
    
    # Standardize rows vectors
    all_rows = _all_rows(df)
    rows_use = _which_rows(rows_use, df)
    
    ### Make dataframe edits
    # Adjustments
    if color_adjustment!=None or color_adj_fxn!=None:
        new_color_by = color_by + "-adj"
        df[new_color_by] = _col(color_by, df, color_adjustment, color_adj_fxn)
        color_by = new_color_by
    # Relabels/reorders
    # if rename_color_groups!=None:
    #     df[color_by] = _rename_and_or_reorder(
    #         df[color_by], reorder = None, relabels = rename_color_groups)
    # if rename_shape_groups!=None:
    #     df[shape_by] = _rename_and_or_reorder(
    #         df[shape_by], reorder = None, relabels = rename_shape_groups)
    # Trim by rows.use, then order if wanted
    Target_data = df.loc[ rows_use ]
    Others_data = df.loc[ [i for i in all_rows if i not in rows_use] ]
    # Reorder for plotting if wanted
    if plot_order != "unordered":
        new_order = order(list(Target_data[color_by]), return_indexes=True)
        if plot_order=="decreasing":
            new_order = reversed(new_order)
        Target_data = Target_data.loc[new_order]
    
    ### Start plot, with data and theming
    fig = (ggplot() +
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
        'data': Target_data,
        'size': size,
        'alpha': opacity}
    
    if color_by!='':
        aes_args['color'] = color_by
        
        if _is_continuous(df[color_by]):
            scale_args= {
                'name': legend_color_title,
                'low': min_color,
                'high': max_color,
                'limits': (
                    [min_value,min(df[color_by])][min_value==None]
                    [max_value,max(df[color_by])][max_value==None])
            }
            if legend_color_breaks!="make":
                scale_args['breaks'] = legend_color_breaks,
            if legend_color_breaks_labels!="make":
                scale_args['labels'] = legend_color_breaks_labels
            fig += scale_colour_gradient(**scale_args)
        else:
            fig += scale_colour_manual(
                name = legend_color_title,
                values = color_panel)
            fig += guides(color = guide_legend(override_aes = {'size':legend_color_size}))
    
    if shape_by!='':
        aes_args['shape'] = shape_by
        fig += scale_shape_manual(
                values = shape_panel,
                name = legend_shape_title)
        fig += guides(shape = guide_legend(override_aes = {'size':legend_shape_size}))
    else:
        geom_args['shape'] = shape_panel[0]
    
    ### Add Data
    # Others_data
    if show_others:
        if (split_show_all_others and split_by!=[]):
            Others_data = _rep_all_data_per_facet(Target_data, Others_data, split_by)
        if Others_data.shape[0] > 0:
            fig += geom_point(
                data = Others_data,
                mapping = aes(x = x_by, y = y_by),
                size = size,
                color = "#E5E5E5")
    # Target_data
    geom_args['mapping'] = aes(**aes_args)
    fig += geom_point(**geom_args)
    
    ### Extra tweaks
    # Faceting
    fig = _add_splitting(fig, split_by, split_nrow, split_ncol, split_adjust)
    # Contours
    if do_contour:
        fig = _add_contours(fig, df, x_by, y_by, contour_color, contour_linetype)
    # Grid lines
    if not show_grid_lines:
        fig += theme(panel_grid_major = element_blank(), panel_grid_minor = element_blank())
    # Legend
    if not legend_show:
        fig = _remove_legend(fig)
    
    if data_out:
        return {
            'plot': fig,
            'Target_data': Target_data,
            'Others_data': Others_data
        }
    return fig
    
def _rep_all_data_per_facet(Target_data, Others_data, split_by):
    '''
    To power showing all data points across all plot facets,
    this function creates a copy of all data points for every facet that will be created.
    It is intended to be used with the "background" 'Others_data' of scatter_plotnine.
    
    'Target_data' a pandas DataFrame of the points which will be the main target of the plot.
        This data frame determines which facets will be shown, so how many replicates of all the data that will be output.
    'Others_data' a pandas DataFrame of the points which were originally not targeted by the upstream function's 'rows_use'.
        Represents the only points that would have been shown in gray without faceting.
        Both 'Target_data' and 'Others_data' should have the same columns as they are splits of the same upstream data frame.
    'split_by' a list of column names (within Target/Others_data) used for faceting.
    
    Details: Target_data and Others_data are recombined into the data frame of all points.
    Then the facets needed are determined.  Others_data in then turned into a bunch of
    row-combined replicates of "all"_data in which each replicate has a given facets' information
    for its 'split_by' columns.  Thus, these points would show up across all facets of the ultimate plot. 
    '''
    
    if split_by==[]:
        return Others_data
        
    all_data = Target_data.append(Others_data).copy()
    
    if len(split_by)==1:
        facets = all_data[split_by[0]]
    if len(split_by)==2:
        facets = [str(a)+str(b) for a,b in zip(all_data[split_by[0]], all_data[split_by[1]])]
    
    Others_data = all_data.loc[[]]
    
    for this_facet in np.unique(facets,True)[1]:
        new_data = all_data.copy()
        # Replace facet info
        for by in split_by:
            new_data[by] = new_data[by][this_facet]
        Others_data = Others_data.append(new_data)
    
    return Others_data

