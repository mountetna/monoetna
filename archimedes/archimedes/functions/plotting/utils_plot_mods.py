from .utils import _is_discrete, _is_continuous, _is_logical, _is_integer, _scale, _all_rows, _which_rows
from .utils_col_getters import _isCol, _col, _colLevels
from plotnine import facet_grid, facet_wrap, theme, aes, geom_density_2d

def _add_splitting(fig, split_by, nrow=None, ncol=None, split_args={}):
    '''
    For plotnine plotting system only
    Adds faceting to go with 'split_by' utilization.
    '''
    if len(split_by)==1:
        if nrow!=None:
            split_args['nrow']=nrow
        if ncol!=None:
            split_args['ncol']=ncol
        fig += facet_wrap(split_by, **split_args)
    if len(split_by)==2:
        fig += facet_grid(split_by, **split_args)
    return fig

def _add_contours(
    fig, data, x_by, y_by, color, linetype = 1):
    '''
    For plotnine plotting system only
    Add contours based on the density of data points in a scatter plot
    '''
    
    return fig + geom_density_2d(
        data = data,
        mapping = aes(x = x_by, y = y_by),
        color = color,
        linetype = linetype,
        na_rm = True)

def _remove_legend(fig):
    '''
    For plotnine plotting system only
    Simply a shorthand for legend removal.
    '''
    return fig + theme(legend_position = "none")
