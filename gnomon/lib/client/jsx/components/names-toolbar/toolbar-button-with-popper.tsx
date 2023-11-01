import React, { useRef } from "react";
import { makeStyles, useTheme, createTheme, ThemeProvider } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import Tooltip from '@material-ui/core/Tooltip';
import Color from 'color';



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


const ToolbarButtonWithPopper = ({ text, iconComponent, variant, color = "default", disabled = false, popperComponent, popperId, onClickOrPopperChange, popperOpen = false }: {
    text: string,
    iconComponent: JSX.Element,
    variant: "full" | "compact",
    color: string,
    disabled?: boolean,
    popperComponent?: JSX.Element,
    popperId?: string,
    onClickOrPopperChange?: (open: boolean) => void,
    popperOpen?: boolean,
}) => {

    const classes = useStyles()
    const anchorEl = useRef(null)
    const baseTheme = useTheme()
    let overrideTheme = baseTheme

    const customColor = ["inherit", "default", "primary", "secondary"].indexOf(color) == -1

    if (customColor) {
        const primary = Color(color)

        overrideTheme = createTheme(baseTheme, {
            palette: {
                primary: {
                    main: primary.string()
                }
            },
            overrides: {
                MuiButton: {
                    containedPrimary: {
                        backgroundColor: primary.alpha(0.3).string(),
                        color: primary.darken(0.5).string(),
                        '&:hover': {
                            backgroundColor: primary.alpha(0.5).string()
                        },
                        '&:disabled': {
                            backgroundColor: primary.alpha(0.15).string()
                        },
                    },
                },
            }
        })
    }

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
            // @ts-ignore
            color={customColor ? "primary" : color}
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
            <ThemeProvider theme={overrideTheme}>
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
            </ThemeProvider>
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