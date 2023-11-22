import React from 'react';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Popper from '@material-ui/core/Popper';
import Grow from '@material-ui/core/Grow';
import ClickAwayListener from '@material-ui/core/ClickAwayListener';
import Paper from '@material-ui/core/Paper';



const useStyles = makeStyles((theme) => ({
    paper: {
        borderRadius: '6px',
        padding: '1em',
    },
    popper: {
        maxWidth: 'min-content',
        '& > :not(:last-child)': {
            marginBottom: '1em',
        }
    },
    buttonsContainer: {
        display: 'flex',
        justifyContent: 'center',
        '& button:not(:last-child)': {
            marginRight: '1em',
        },
    },
}));


const ConfirmationPopper = ({ text, open, onConfirm, onClose, className, anchorRef }: {
    text?: string
    open: boolean,
    onConfirm: (confirmed: boolean) => void,
    onClose: () => void,
    className?: string
    anchorRef: React.RefObject<HTMLElement>,
}) => {
    const classes = useStyles();

    return (
        <Popper
            open={open}
            anchorEl={anchorRef.current}
            placement='bottom'
            role={undefined}
            transition
        >
            {({ TransitionProps, placement }) => (
                <Grow
                    {...TransitionProps}
                    style={{ transformOrigin: placement === 'bottom' ? 'center top' : 'center bottom' }}
                    exit={false}
                >
                    <Paper variant="outlined" className={classes.paper}>
                        <ClickAwayListener onClickAway={onClose}>

                            <div className={`${classes.popper} ${className != undefined ? className : ''}`}>
                                {text?.length && <div>{text}</div>}

                                <div className={classes.buttonsContainer}>
                                    <Button
                                        variant="contained"
                                        color="secondary"
                                        disableRipple
                                        disableElevation
                                        aria-label='cancel-confirmation'
                                        onClick={() => onConfirm(false)}
                                    >
                                        Cancel
                                    </Button>

                                    <Button
                                        variant="contained"
                                        color="primary"
                                        disableRipple
                                        disableElevation
                                        aria-label='confirm-confirmation'
                                        onClick={() => onConfirm(true)}
                                    >
                                        Confirm
                                    </Button>
                                </div>
                            </div>

                        </ClickAwayListener>
                    </Paper>
                </Grow>
            )}
        </Popper>
    );
};


export default ConfirmationPopper;