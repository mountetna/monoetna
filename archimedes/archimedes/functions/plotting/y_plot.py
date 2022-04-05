import plotly.express as px

from plotnine import ggplot, aes, theme_bw, ggtitle, geom_jitter, geom_violin, geom_boxplot, scale_fill_manual, xlab, ylab, facet_wrap, facet_grid, position_jitterdodge, position_dodge, geom_hline, element_text, theme, scale_y_continuous, coord_cartesian

from .utils import _leave_default_or_none, _which_rows
from .colors import colors

def y_plotly(
    data_frame,
    x_by,
    y_by,
    color_by = "make",
    plots = ["violin", "box"],
    px_args: dict = {},
    rows_use = None,
    y_scale = "linear",
    color_panel: list = colors,
    x_title="make",
    y_label="make",
    plot_title="make",
    legend_title ="make"):
    """
    Produces violin and/or box plots using plotly.express.violin, or plotly.express.box, based on the pandas 'data_frame' given.
    'x_by', 'y_by' should indicate columns of 'data_frame' to use for x/y axes data. 'x_by' should denote discrete groupings, while 'y-by' should denote numerical data.
    'color_by' can indicate a column of 'data_frame' to use for establishing fill colors of violin/boxplots and can be used to create subgroupings, or highlight supergroupings, of 'x_by' groups. By default, if left as "make", it will be passed the 'x_by' data.
    'plots' sets which data representations to use, and should be a list containing "violin", "box", or both of these.
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.violin' or 'plotly.express.box' call.
    'color_panel' (string list) sets the colors to use for violin/boxplot fills.
    'plot_title', 'legend_title', 'x_title', and 'y_title' set titles.
    'rows_use'
    'y_scale', String, 'linear', 'log10', or 'log10(val+1)'. Controls whether this axes should be log scaled, and if so, whether 1 should be added to all values first in order to let zeros be okay to plot.
    """
    
    # Parse dependent defaults
    color_by = _leave_default_or_none(color_by, x_by)
    x_title = _leave_default_or_none(x_title, x_by)
    y_label = _leave_default_or_none(y_label, y_by)
    plot_title = _leave_default_or_none(plot_title, y_label)
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
        xaxis_title=x_title,
        yaxis_title=y_label,
        legend_title_text=legend_title,
        legend= {'itemsizing': 'constant'}
    )
    return fig

def y_plotnine(
    data_frame,
    x_by,
    y_by,
    color_by = "make",
    plots = ["violin", "box", "jitter"],
    split_by = [],
    # indexes_use = None,
    color_panel: list = colors,
    boxplot_fill = True,
    boxplot_width = 0.2,
    boxplot_color = "black",
    boxplot_show_outliers = "make",
    boxplot_position_dodge = "make",
    jitter_width = 0.2,
    jitter_size = 1,
    jitter_color = "black",
    jitter_position_dodge = "make",
    violin_lineweight = 1,
    violin_width = 1,
    violin_scale = "area",
    violin_quantiles = None,
    add_line = None,
    line_linetype = "dashed",
    line_color = "black",
    y_scale = scale_y_continuous,
    y_breaks = None,
    y_min = None,
    y_max = None,
    x_labels_rotate = True,
    x_title="make",
    y_title="make",
    plot_title="make",
    legend_title ="make",
    plot_theme = theme_bw()
    ):
    """
    NOT UPDATED
    Produces violin and/or box plots using plotly.express.violin, or plotly.express.box, based on the pandas 'data_frame' given.
    'x_by', 'y_by' should indicate columns of 'data_frame' to use for x/y axes data. 'x_by' should denote discrete groupings, while 'y-by' should denote numerical data.
    'color_by' can indicate a column of 'data_frame' to use for establishing fill colors of violin/boxplots and can be used to create subgroupings, or highlight supergroupings, of 'x_by' groups. By default, if left as "make", it will be passed the 'x_by' data.
    'plots' sets which data representations to use, and should be a list containing "violin", "box", or both of these.
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.violin' or 'plotly.express.box' call.
    'color_panel' (string list) sets the colors to use for violin/boxplot fills.
    'plot_title', 'legend_title', 'x_title', and 'y_title' set titles.
    """
    
    # Parse dependent defaults
    color_by = _leave_default_or_none(color_by, x_by)
    x_title = _leave_default_or_none(x_title, x_by)
    y_title = _leave_default_or_none(y_title, y_by)
    plot_title = _leave_default_or_none(plot_title, None)
    legend_title = _leave_default_or_none(legend_title, color_by)
    boxplot_show_outliers = _leave_default_or_none(boxplot_show_outliers, "jitter" not in plots)
    boxplot_position_dodge = _leave_default_or_none(boxplot_position_dodge, violin_width)
    jitter_position_dodge = _leave_default_or_none(jitter_position_dodge, boxplot_position_dodge)    
    
    # # Trim df if requested
    # if indexes_use!=None:
    #     indexes_use=_which_indexes(indexes_use, data_frame)
    
    ### Start Plot, with data and theming
    fig = ( ggplot(data_frame, aes(x = x_by, y = y_by, fill = color_by)) +
        plot_theme +
        scale_fill_manual(values = color_panel, name = legend_title) +
        xlab(x_title) + ylab(y_title)
        )
    
    if plot_title!=None:
        fig += ggtitle(plot_title)
        
    if y_breaks!=None:
        fig += y_scale(breaks = y_breaks)
    else:
        fig += y_scale()
        
    if y_min==None:
        y_min = min(data_frame[y_by])
    if y_max==None:
        y_max = max(data_frame[y_by])
    fig += coord_cartesian(ylim=(y_min,y_max))
    
    # Add data
    for ptype in plots:
        if ptype == "violin":
            violin_args = {
                'mapping': aes(fill=color_by),
                'width': violin_width,
                'scale': violin_scale,
                'weight': violin_lineweight,
                'na_rm': True
            }
            if violin_quantiles!=None:
                violin_args['draw_quantiles'] = violin_quantiles
            fig += geom_violin(**violin_args)
        if ptype == "box":
            boxplot_args = {
                'mapping': aes(fill=color_by),
                'width': boxplot_width,
                'color': boxplot_color,
                'alpha': [0,1][boxplot_fill],
                'outlier_shape': ['','o'][boxplot_show_outliers],
                'position': position_dodge(width = boxplot_position_dodge),
                'na_rm': True
            }
            fig += geom_boxplot(**boxplot_args)
        if ptype == "jitter":
            fig += geom_jitter(
                position = position_jitterdodge(
                      jitter_width = jitter_width,
                      jitter_height = 0,
                      dodge_width = jitter_position_dodge,
                ),
                shape = '.',
                size = jitter_size,
                color = jitter_color,
                na_rm = True
            )
    
    ### Extra tweaks
    # Faceting
    if len(split_by)==1:
        fig += facet_wrap(split_by)
    if len(split_by)==2:
        fig += facet_grid(split_by)
    # Horizontal Lines
    if add_line!=None:
        fig += geom_hline(yintercept=add_line, linetype= line_linetype, color = line_color)
    # Rotated labels
    if x_labels_rotate:
        fig += theme(axis_text_x= element_text(angle=45, hjust = 1, vjust = 1))
    
    return fig

