import React, { useState, useRef } from "react";
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import AddCircleOutlineIcon from "@material-ui/icons/AddCircleOutline";


const useStyles = makeStyles((theme) => ({
    addNameContainer: {
        display: "inline-block",
    },
}));


const AddNamesButton = ({ rules, clickAddHandler }: { rules: string[], clickAddHandler: (rule: string) => any }) => {
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
        clickAddHandler(rule);
        handleClose();
    };

    return (
        <div className={classes.addNameContainer}>
            <Button
                startIcon={<AddCircleOutlineIcon />}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label="Add Names"
                aria-haspopup="true"
                aria-controls={open ? "add-names-rules-menu" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
            >
                Add Names
            </Button>
            <Popper
                open={open}
                anchorEl={anchorEl.current}
                placement='bottom'
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
                                <MenuList autoFocusItem={open} id="add-names-rules-menu">
                                    {
                                        rules.map((rule) =>
                                            <MenuItem
                                                onClick={() => handleClickRule(rule)}
                                                key={rule}
                                                disableRipple
                                            >
                                                {rule}
                                            </MenuItem>
                                        )
                                    }
                                </MenuList>
                            </ClickAwayListener>
                        </Paper>
                    </Grow>
                )}
            </Popper>
        </div>
    )
};

export default AddNamesButton;