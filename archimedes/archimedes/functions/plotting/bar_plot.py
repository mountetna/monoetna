import plotly.express as px
import pandas as pd

from .utils import _leave_default_or_none, _which_rows
from .colors import colors
from ..list import unique, order

def bar_plotly(
    data_frame,
    x_by,
    y_by,
    scale_by = 'fraction',
    px_args: dict = {},
    rows_use = None,
    color_panel: list = colors,
    xlab = "make",
    ylab = "make",
    plot_title = "make",
    legend_title = "make"):
    """
    Produces a stacked bar plot using 'plotly.express.bar' based on the pandas 'data_frame' given, which describes the composition of 'y_by' sets within 'x_by' groupings.
    'x_by', 'y_by' should indicate columns of 'data_frame' to use for x/y axes data.  Both must point to discrete data.
    'scale_by' can be either "fraction" or "counts" and sets whether the y-axis of the plot will ultimately reflect proportions of 'y_by'-value make-up within each 'x_by' group versus raw counts of 'y_by'-values per group. 
    'px_args' should be a dictionary of additional bits to send in the 'plotly.express.bar' call.
    'color_panel' (string list) sets the colors to assign to the distinct 'y_by'-values.
    'plot_title', 'legend_title', 'xlab', and 'ylab' set titles.
    
    Some additional details: Prior to making a plot, the composition of the 'y_by' values associated with each 'x_by' grouping is calculated and scaled based on the method indicated by the 'scale_by' parameter.
    The resulting summary dataframe is then passed to plotly for plot creation.
    """

    # Parse dependent defaults
    xlab = _leave_default_or_none(xlab, x_by)
    ylab = _leave_default_or_none(ylab, y_by + " " + scale_by) # A little different for this function
    plot_title = _leave_default_or_none(plot_title, ylab + " per " + x_by)
    legend_title = _leave_default_or_none(legend_title, y_by) # A little different for this function
    
    # data_frame edits
    df = data_frame.copy()
    rows_use = _which_rows(rows_use, df)
    df = df.loc[rows_use]

    ### Generate a composition summary dataframe from the input df.
    summary_df = pd.DataFrame(df[[x_by, y_by]].value_counts())
    summary_df = summary_df.reset_index().rename(columns={x_by: 'xgroups', y_by: 'yvals', 0 : "counts"})
    fractions = None
    for k in summary_df.xgroups.unique(): # for each unique value of the xgroup column, refered to now as k
        counts = summary_df[summary_df.xgroups == k].loc[:,'counts'] # grabbing all the counts values for the current xgroup value level
        this = counts/sum(counts) # calculate the percent that each counts value makes up
        if fractions is None:
            fractions = this
        else:
            fractions = fractions.append(this)
    summary_df['fraction'] = fractions
    
    # Add to px_args and convert from our variables names to px.violin/bar variables names.
    px_args['data_frame'] = summary_df
    px_args['x'] = "xgroups"
    px_args['y'] = scale_by
    px_args['color'] = "yvals"
    px_args['color_discrete_sequence'] = color_panel
    
    ### Make plot
    fig = px.bar(**px_args)

    # Tweaks
    fig.update_layout(
        title_text=plot_title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        legend_title_text=legend_title,
        legend= {'itemsizing': 'constant'}
    )

    return fig

