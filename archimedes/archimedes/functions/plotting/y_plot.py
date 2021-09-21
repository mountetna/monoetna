import plotly.express as px

from .utils import _default_to_if_make_and_logic
from .colors import colors

def y_plotly(
    data_frame,
    x_by,
    y_by,
    color_by = "make",
    plots = ["violin", "box"],
    xlab="make",
    ylab="make",
    plot_title="make",
    legend_title ="make"):
    
    # Parse dependent defaults
    color_by = _default_to_if_make_and_logic(color_by, x_by)
    xlab = _default_to_if_make_and_logic(xlab, x_by)
    ylab = _default_to_if_make_and_logic(ylab, y_by)
    plot_title = _default_to_if_make_and_logic(plot_title, ylab)
    legend_title = _default_to_if_make_and_logic(legend_title, color_by)
    
    # Convert from our variables names to px.violin/bar variables names. 
    the_atributes = {
        "data_frame": data_frame,
        "x": x_by,
        "y": y_by,
        "color": color_by,
        "color_discrete_sequence": colors
        }
    
    # Make Plot
    if "violin" in plots:
        #use conditional here for checking if box is in plots
        the_atributes["box"]= "box" in plots
        fig = px.violin(
            **the_atributes
            )
    elif "box" in plots:
        fig = px.box(
            **the_atributes
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

