import plotly.express as px
import pandas as pd

from .utils import _default_to_if_make_and_logic
from .colors import colors
from ..list import unique, order

def bar_plotly(
    data_frame,
    x_by,
    y_by,
    scale_by = 'counts',
    xlab = "make",
    y_lab = "make",
    main = "make",
    legend_title = "make"):

    # Parse dependent defaults
    xlab = _default_to_if_make_and_logic(xlab, x_by)
    ylab = _default_to_if_make_and_logic(ylab, y_by + " " + scale_by)
    plot_title = _default_to_if_make_and_logic(plot_title, ylab)
    legend_title = _default_to_if_make_and_logic(legend_title, y_by + "\n" + scale_by)

    ### Generate a composition summary dataframe from the input df.
    summary_df = pd.DataFrame(data_frame[[x_by, y_by]].value_counts())
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
    ### Make plot
    fig = px.bar(
        data_frame = summary_df,
        x = "xgroups",
        y = scale_by,
        color="yvals",
        color_discrete_sequence = colors
    )

    # Tweaks
    fig.update_layout(
        title_text=plot_title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        legend= {'itemsizing': 'constant'}
    )
    fig.update_coloraxes(colorbar_title_text=legend_title)

    return fig

