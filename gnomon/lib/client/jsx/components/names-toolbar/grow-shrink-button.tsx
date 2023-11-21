import React from 'react';
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';



const GrowShrinkButton = ({ small, onClick, className }: {
    small: boolean,
    onClick: () => void,
    className?: string
}) => {
    return (
        <ToolbarButtonWithPopper
            text={small ? 'Grow' : 'Shrink'}
            iconComponent={small ? <UnfoldMoreOutlinedIcon /> : <UnfoldLessOutlinedIcon />}
            variant={small ? 'compact' : 'full'}
            color="secondary"
            onClickOrPopperChange={onClick}
            className={className}
        />
    );
};

export default GrowShrinkButton;