import React, { useState, useEffect, ReactNode, useRef } from "react";
import Button from "@material-ui/core/Button";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown';

import { CreateName, Rule, Token, TokenValue } from "../../models";




const TokenMenu = ({ token, initialValue }: { token: Token, initialValue?: TokenValue }) => {
    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null)
    const [value, setValue] = useState<TokenValue>(initialValue)

    const handleToggle = () => {
        setOpen(prev => !prev);
    };
    const handleClose = () => {
        setOpen(false);
    };
    const handleClickTokenValue = (value: TokenValue) => {
        console.log(`setting value ${value} for token ${token}`)
        setValue(value)
        handleClose();
    };

    return (
        <div className="select-token-value-container">
            <Button
                startIcon={KeyboardArrowDownIcon}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label={`Select Token Value for "${token.label}"`}
                aria-haspopup="true"
                // aria-controls={open ? "add-names-rules-menu" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
            >
                {value ? value.name : token.name}
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
                                <MenuList autoFocusItem={open} id="add-names-rules-menu">
                                    <MenuItem key={token.label}>
                                        Choose a {token.label}
                                    </MenuItem>
                                    {
                                        token.values.map((val) =>
                                            <MenuItem
                                                onClick={() => handleClickTokenValue(val)}
                                                key={val}
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



const NameElementsEditor = ({ name, rule }: { name: CreateName, rule: Rule }) => {
    return (
        <span className="name-elements-editor">
            {
                rule.tokens.map((el, idx) => {
                    // it's a token
                    if (el.hasOwnProperty("name")) {
                        return <span key={idx}></span>
                    }
                })
            }
        </span>
    )
}



const NameComposer = ({ name, rule }: { name: CreateName, rule: Rule }) => {
    return (
        <div className="name-composer-container">
            <span className="name-composer-tools"></span>
            <NameElementsEditor name={name} rule={rule} />
        </div>
    )
}


export default NameComposer;