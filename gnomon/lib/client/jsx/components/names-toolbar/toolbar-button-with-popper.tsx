import React, { useRef } from 'react';
import { makeStyles, useTheme, ThemeProvider } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Popper from '@material-ui/core/Popper';
import Grow from '@material-ui/core/Grow';
import ClickAwayListener from '@material-ui/core/ClickAwayListener';
import Paper from '@material-ui/core/Paper';
import Tooltip from '@material-ui/core/Tooltip';
import { createCustomColorButtonTheme } from '../../utils/theme';



const useStyles = makeStyles((theme) => ({
    container: {
        display: 'inline-block',
    },
    button: {
        height: '4em',
        lineHeight: '1.5em',
    },
    popperContainer: {
        borderRadius: '6px',
    },
}));


const ToolbarButtonWithPopper = ({ text, iconComponent, variant, color = 'default', disabled = false, popperComponent, popperId, onClickOrPopperChange, popperOpen = false, className, buttonRef }: {
    text: string,
    iconComponent: JSX.Element,
    variant: 'full' | 'compact',
    color: string,
    disabled?: boolean,
    popperComponent?: JSX.Element,
    popperId?: string,
    onClickOrPopperChange?: (open: boolean) => void,
    popperOpen?: boolean,
    className?: string,
    buttonRef?: React.MutableRefObject<null>
}) => {

    const _anchorRef = useRef(null);
    const anchorRef = buttonRef || _anchorRef;

    const classes = useStyles();
    const baseTheme = useTheme();
    const customColor = ['inherit', 'default', 'primary', 'secondary'].indexOf(color) == -1;
    const theme = customColor ? createCustomColorButtonTheme(color, baseTheme) : baseTheme;

    const handleToggle = () => {
        if (onClickOrPopperChange != undefined) {
            onClickOrPopperChange(!popperOpen);
        }
    };

    const handleClose = () => {
        if (onClickOrPopperChange != undefined) {
            onClickOrPopperChange(false);
        }
    };

    const renderButton = () => {
        return <Button
            startIcon={variant == 'full' ? iconComponent : undefined}
            onClick={handleToggle}
            ref={anchorRef}
            // @ts-ignore
            color={customColor ? 'primary' : color}
            aria-label={text}
            aria-haspopup={popperComponent != undefined ? 'true' : 'false'}
            aria-controls={popperOpen ? popperId : undefined}
            size={'medium'}
            className={classes.button}
            disableRipple
            disableFocusRipple
            disableElevation
            disabled={disabled}
        >
            {variant == 'full' ? text : iconComponent}
        </Button>;
    };

    console.log(className);

    return (
        <div className={`${classes.container} ${className != undefined ? className : ''}`}>
            <ThemeProvider theme={theme}>
                {
                    // https://v4.mui.com/components/tooltips/#disabled-elements
                    // using <span> to wrap makes <Tooltip> positioning inconsistent
                    variant == 'full' || disabled == true
                        ? renderButton()
                        : (
                            <Tooltip title={text} placement="top">
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
                        anchorEl={anchorRef.current}
                        placement='bottom'
                        role={undefined}
                        transition
                    >
                        {({ TransitionProps, placement }) => (
                            <Grow
                                {...TransitionProps}
                                style={{ transformOrigin: placement === 'bottom' ? 'center top' : 'center bottom' }}
                            >
                                <Paper variant="outlined" className={classes.popperContainer}>
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
    );
};

export default ToolbarButtonWithPopper;