import React, { useState, useRef } from "react";
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import Tooltip from '@material-ui/core/Tooltip';



const useStyles = makeStyles((theme) => ({
    container: {
        display: "inline-block",
        maxWidth: "8em",
        marginRight: "1.25em",
        "&:last-child": {
            marginRight: "0",
        },
    },
    button: {
        height: "4em",
        lineHeight: "1.5em",
    },
}));


const ToolbarButtonWithPopper = ({ text, iconComponent, variant, color, disabled = false, popperComponent, popperId, onClickOrPopperChange, popperOpen = false }: {
    text: string,
    iconComponent: JSX.Element,
    variant: "full" | "compact",
    color?: "inherit" | "default" | "primary" | "secondary",
    disabled?: boolean,
    popperComponent?: JSX.Element,
    popperId?: string,
    onClickOrPopperChange?: (open: boolean) => void,
    popperOpen?: boolean,
}) => {

    const classes = useStyles()
    const anchorEl = useRef(null)

    const handleToggle = () => {
        if (onClickOrPopperChange != undefined) {
            onClickOrPopperChange(!popperOpen)
        }
    };

    const handleClose = () => {
        if (onClickOrPopperChange != undefined) {
            onClickOrPopperChange(false);
        }
    };

    const renderButton = () => {
        return <Button
            startIcon={variant == "full" ? iconComponent : undefined}
            onClick={handleToggle}
            ref={anchorEl}
            color={color}
            aria-label={text}
            aria-haspopup={popperComponent != undefined ? "true" : "false"}
            aria-controls={popperOpen ? popperId : undefined}
            size={variant == "full" ? "medium" : "small"}
            className={classes.button}
            disableRipple
            disableFocusRipple
            disableElevation
            disabled={disabled}
        >
            {variant == "full" ? text : iconComponent}
        </Button>
    }

    return (
        <div className={classes.container}>
            {
                variant == "full" || disabled == true
                    ? renderButton()
                    : (
                        <Tooltip title={text} placement="top">
                            {/* https://v4.mui.com/components/tooltips/#disabled-elements */}
                            {renderButton()}
                        </Tooltip>
                    )
            }
            {
                popperComponent != undefined
                && (
                    <Popper
                        open={popperOpen}
                        anchorEl={anchorEl.current}
                        placement='bottom'
                        role={undefined}
                        transition
                    >
                        {({ TransitionProps }) => (
                            <Grow
                                {...TransitionProps}
                                style={{ transformOrigin: "center top" }}
                            >
                                <Paper variant="outlined">
                                    <ClickAwayListener onClickAway={handleClose}>
                                        {popperComponent}
                                    </ClickAwayListener>
                                </Paper>
                            </Grow>
                        )}
                    </Popper>
                )
            }
        </div>
    )
};

export default ToolbarButtonWithPopper;