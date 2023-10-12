import React, { useState, useRef } from "react";
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import ButtonBase from "@material-ui/core/ButtonBase";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown';

import { Rule, Token, TokenValue, UNSET_TOKEN_VALUE } from "../../../models";
import { selectTokenValueLocalIdsWithTokenName, selectTokenValuesByLocalId } from "../../../selectors/rules";



const useStyles = makeStyles((theme) => ({
    selectValueContainer: {
        display: "inline-block",
    },
    currentSelectValue: {
        fontWeight: "bold"
    },
    unset: {
        color: "red"
    },
}));


interface SelectValue {
    name: string
    label?: string
}


// TODO: change to actual select
const SelectBase = <T extends SelectValue>({ values, value, label, placeholder, onSetValue }:
    {
        values: T[],
        value?: T,
        label: string,
        placeholder: string,
        onSetValue: (value: T) => void,
    }) => {

    const classes = useStyles()
    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null)

    const canChooseValue = () => {
        if (value == undefined) {
            return true
        }
        if (values.length == 0) {
            return false
        }
        if (values.length == 1) {
            return value == undefined || value != values[0]
        }
        if (values.length > 1) {
            return true
        }
    }

    const handleToggle = () => {
        setOpen(prev => !prev);
    };

    const handleClose = () => {
        setOpen(false);
    };

    const handleClickMenuItem = (menuItemValue: T) => {
        if (menuItemValue != value) {
            onSetValue(menuItemValue)
        }
        handleClose()
    };

    return (
        <div className={classes.selectValueContainer}>
            <ButtonBase
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label={`Select Value for "${label}"`}
                aria-haspopup="true"
                aria-controls={open ? "select-value-container" : undefined}
                disableRipple
                disableTouchRipple
                disabled={!canChooseValue()}
            >
                <span
                    className={classes.currentSelectValue + " " + (!value ? `${classes.unset}` : "")}
                >
                    {canChooseValue() && <KeyboardArrowDownIcon />}
                    {value ? value.name : placeholder}
                </span>
            </ButtonBase>
            <Popper
                open={open}
                anchorEl={anchorEl.current}
                placement="bottom"
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
                                <MenuList autoFocusItem={open} id="select-value-container">
                                    <div key={label}>
                                        Choose a {label}
                                    </div>
                                    {
                                        values.map(val =>
                                            <MenuItem
                                                onClick={() => handleClickMenuItem(val)}
                                                key={val.name}
                                                disableRipple
                                            >
                                                {val.name}{val.label ? ` - ${val.label}` : ""}
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


export const RuleSelect = ({ values, value, label, placeholder, onSetRule }:
    {
        values: Rule[],
        value?: Rule,
        label: string,
        placeholder: string,
        onSetRule: (value: Rule) => void,
    }
) => {

    return (
        <SelectBase
            values={values}
            value={value}
            label={label}
            placeholder={placeholder}
            onSetValue={onSetRule}
        />
    )
}


export const TokenSelect = ({ token, value, onSetTokenValue, includeUnsetAsValue = false }:
    { token: Token, value?: TokenValue, onSetTokenValue: (value: TokenValue) => void, includeUnsetAsValue?: boolean }) => {

    const tokenValuesByLocalId: Record<string, TokenValue> = useSelector(selectTokenValuesByLocalId)
    const tokenValues: TokenValue[] = useSelector(selectTokenValueLocalIdsWithTokenName(token.name))
        .map(tvLocalId => tokenValuesByLocalId[tvLocalId])

    if (includeUnsetAsValue === true) {
        tokenValues.unshift(UNSET_TOKEN_VALUE)
    }

    return (
        <SelectBase
            values={tokenValues}
            value={value}
            label={token.label}
            placeholder={token.name}
            onSetValue={onSetTokenValue}
        />
    )
}

