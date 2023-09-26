import React, { useState, useRef } from "react";
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import LibraryAddOutlinedIcon from "@material-ui/icons/LibraryAddOutlined";

import { CreateName, Rule } from "../../models";


const useStyles = makeStyles((theme) => ({
    addFromSelectionContainer: {
        display: "inline-block",
    },
}));


const AddFromSelectionButton = ({ selection, rules, clickAddHandler }:
    {
        selection: CreateName[],
        rules: Rule[],
        clickAddHandler: (names: CreateName[], tokenValue: string, start: number, finish: number) => any
    }) => {

    const classes = useStyles()
    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null)

    const handleToggle = () => {
        setOpen(prev => !prev);
    };
    const handleClose = () => {
        setOpen(false);
    };
    const handleClickRule = (rule: string) => {
        clickAddHandler(selection);
        handleClose();
    };

    return (
        <div className={classes.addFromSelectionContainer}>
            <Button
                startIcon={<LibraryAddOutlinedIcon />}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label="Add from Selection"
                aria-haspopup="true"
                aria-controls={open ? "add-from-selection-dialogue" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
            >
                Add from Selection
            </Button>
            <Popper
                open={open}
                anchorEl={anchorEl.current}
                placement="bottom"
                role={undefined}
                transition
                disablePortal
            >
                {({ TransitionProps }) => (
                    <Grow
                        {...TransitionProps}
                        style={{ transformOrigin: "center top" }}
                    >
                        <Paper variant="outlined">
                            <ClickAwayListener onClickAway={handleClose}>
                                <div>hi</div>
                            </ClickAwayListener>
                        </Paper>
                    </Grow>
                )}
            </Popper>
        </div>
    )
};

export default AddFromSelectionButton;