import React from 'react';
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import { makeStyles } from '@material-ui/core/styles';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';



const useStyles = makeStyles((theme) => ({
    container: {
        '&.iconVertical svg': {},
        '&.iconHorizontal svg': {
            rotate: '90deg',
        },
    },
}));


const GrowShrinkButton = ({ small, onClick, className, growText = 'Grow', shrinkText = 'Shrink', iconDirection = 'vertical' }: {
    small: boolean,
    onClick: () => void,
    className?: string,
    growText?: string,
    shrinkText?: string,
    iconDirection?: 'vertical' | 'horizontal',
}) => {
    const classes = useStyles();

    return (
        <ToolbarButtonWithPopper
            text={small ? growText : shrinkText}
            iconComponent={small ? <UnfoldMoreOutlinedIcon /> : <UnfoldLessOutlinedIcon />}
            variant={small ? 'compact' : 'full'}
            color="secondary"
            onClickOrPopperChange={onClick}
            className={
                classes.container
                + (className ? ` ${className}` : '')
                + (iconDirection == 'horizontal' ? ' iconHorizontal' : 'iconVertical')
            }
        />
    );
};

export default GrowShrinkButton;