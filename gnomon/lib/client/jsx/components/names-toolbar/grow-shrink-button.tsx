import React from "react";
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import ToolbarButtonWithPopper from "./toolbar-button-with-popper";



const GrowShrinkButton = ({ small, onClick }: {
    small: boolean,
    onClick: () => void,
}) => {
    return (
        <ToolbarButtonWithPopper
            text={small ? "Grow" : "Shrink"}
            iconComponent={small ? <UnfoldMoreOutlinedIcon /> : <UnfoldLessOutlinedIcon />}
            variant={small ? "compact" : "full"}
            color="primary"
            onClickOrPopperChange={onClick}
        />
    )
};

export default GrowShrinkButton;