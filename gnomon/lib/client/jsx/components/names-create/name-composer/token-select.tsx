import React, { useState, useRef } from "react";
import { makeStyles } from '@material-ui/core/styles';
import ButtonBase from "@material-ui/core/ButtonBase";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown';

import { Token, TokenValue } from "../../../models";



const useStyles = makeStyles((theme) => ({
    selectTokenValueContainer: {
        display: "inline-block",
    },
    currentTokenValue: {
        fontWeight: "bold"
    },
    unset: {
        color: "red"
    },
}));


// TODO: change to actual select
const TokenSelect = ({ token, value, handleSetTokenValue }:
    { token: Token, value?: TokenValue, handleSetTokenValue: (value: TokenValue) => void }) => {

    const classes = useStyles()
    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null)

    const handleToggle = () => {
        setOpen(prev => !prev);
    };
    const handleClose = () => {
        setOpen(false);
    };
    const handleClickTokenValue = (tokenValue: TokenValue) => {
        if (tokenValue != value) {
            handleSetTokenValue(tokenValue)
        }
        handleClose()
    };

    return (
        <div className={classes.selectTokenValueContainer}>
            <ButtonBase
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label={`Select Token Value for "${token.label}"`}
                aria-haspopup="true"
                aria-controls={open ? "select-token-value-container" : undefined}
                disableRipple
                disableTouchRipple
            >
                <span
                    className={classes.currentTokenValue + " " + (!value ? `${classes.unset}` : "")}
                >
                    <KeyboardArrowDownIcon />
                    {value ? value.name : token.name}
                </span>
            </ButtonBase>
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
                                <MenuList autoFocusItem={open} id="select-token-value-container">
                                    <div key={token.label}>
                                        Choose a {token.label}
                                    </div>
                                    {
                                        token.values.map((val) =>
                                            <MenuItem
                                                onClick={() => handleClickTokenValue(val)}
                                                key={val.name}
                                                disableRipple
                                            >
                                                {val.name} - {val.label}
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
}


export default TokenSelect