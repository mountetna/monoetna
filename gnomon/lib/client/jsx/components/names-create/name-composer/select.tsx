import React, { useState, useRef } from 'react';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import MenuList from '@material-ui/core/MenuList';
import MenuItem from '@material-ui/core/MenuItem';
import Popper from '@material-ui/core/Popper';
import Grow from '@material-ui/core/Grow';
import ClickAwayListener from '@material-ui/core/ClickAwayListener';
import Paper from '@material-ui/core/Paper';
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown';

import { Rule, Token, TokenValue, UNSET_TOKEN_VALUE } from '../../../models';
import { selectTokenValueLocalIdsWithTokenName, selectTokenValuesByLocalId } from '../../../selectors/rules';
import { useAppSelector as useSelector } from '../../../hooks';



const useStyles = makeStyles((theme) => ({
    selectValueContainer: {
        display: 'inline-flex',
        '&:not(:last-child)': {
            marginRight: '0.4em',
        },
    },
    currentSelectValue: {
        fontWeight: 'bold',
        background: 'none',
        padding: '0',
        minWidth: 'unset',
        transition: 'background 0s',
        '&:disabled, &:hover, &:active': {
            background: 'none',
            color: 'unset',
        },
        '&.highlighted': {
            background: 'yellow',
        },
        '&:not(:disabled)': {
            // account for dropdown svg icon
            marginLeft: '-0.2em',
        },
        '& .MuiButton-label': {
            lineHeight: '1em',
            fontSize: '16px',
        },
        '& .MuiButton-startIcon': {
            margin: '0',
            opacity: '0.25',
            transition: 'opacity 0.2s ease-in',
        },
        '&:hover .MuiButton-startIcon, &:active .MuiButton-startIcon': {
            opacity: '1',
        },
    },
    unset: {
        color: 'red !important',
        '& .MuiButton-startIcon': {
            opacity: '1',
        },
    },
    nullSelectValue: {
        fontStyle: 'italic',
    },
}));


interface SelectValue {
    name: string
    label?: string
}


// TODO: change to actual select
const SelectBase = <T extends SelectValue>({ values, value, label, placeholder, className, includeNullValue = false, highlight, onSetValue }: {
    values: T[],
    value?: T,
    label: string,
    placeholder: string,
    className?: string,
    includeNullValue?: boolean,
    highlight: boolean,
    onSetValue: (value?: T) => void,
}) => {

    const classes = useStyles({ highlight });

    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null);

    const canChooseValue = () => {
        if (value == undefined) {
            return true;
        }
        if (values.length == 0) {
            return false;
        }
        if (values.length == 1) {
            return value == undefined || value != values[0];
        }
        if (values.length > 1) {
            return true;
        }
        return false;
    };

    const handleToggle = () => {
        setOpen(prev => !prev);
    };

    const handleClose = () => {
        setOpen(false);
    };

    const handleClickMenuItem = (menuItemValue?: T) => {
        if (menuItemValue != value) {
            onSetValue(menuItemValue);
        }
        handleClose();
    };

    return (
        <div className={`${classes.selectValueContainer} ${className != undefined ? className : ''}`}>
            <Button
                startIcon={canChooseValue() ? <KeyboardArrowDownIcon /> : undefined}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label={`Select Value for ${label}`}
                aria-haspopup="true"
                className={
                    classes.currentSelectValue + (!value ? ` ${classes.unset}` : '') + (highlight ? ' highlighted' : '')
                }
                disableRipple
                disableTouchRipple
                disabled={!canChooseValue()}
            >
                {value ? value.name : placeholder}
            </Button>
            <Popper
                open={open}
                anchorEl={anchorEl.current}
                placement="bottom"
                role={undefined}
                transition
            >
                {({ TransitionProps, placement }) => (
                    <Grow
                        {...TransitionProps}
                        style={{ transformOrigin: placement === 'bottom' ? 'center top' : 'center bottom' }}
                    >
                        <Paper variant="outlined">
                            <ClickAwayListener onClickAway={handleClose}>
                                <MenuList autoFocusItem={open}>
                                    <MenuItem key={label} disabled>
                                        <em>Choose a {label}</em>
                                    </MenuItem>
                                    {includeNullValue && <MenuItem
                                        onClick={() => handleClickMenuItem()}
                                        key="none"
                                        disableRipple
                                        className={classes.nullSelectValue}
                                    >
                                        none
                                    </MenuItem>}
                                    {
                                        values.map(val =>
                                            <MenuItem
                                                onClick={() => handleClickMenuItem(val)}
                                                key={val.name}
                                                disableRipple
                                            >
                                                {val.name + (val.label ? ` - ${val.label}` : '')}
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
    );
};


export const RuleSelect = ({ values, value, label, placeholder, className, onSetRule }:
    {
        values: Rule[],
        value?: Rule,
        label: string,
        placeholder: string,
        className?: string,
        onSetRule: (value?: Rule) => void,
    }
) => {

    return (
        <SelectBase
            values={values}
            value={value}
            label={label}
            placeholder={placeholder}
            className={className}
            onSetValue={onSetRule}
            highlight={false}
        />
    );
};


export const TokenSelect = ({ token, value, className, onSetTokenValue, includeUnsetAsValue = false, highlight }: {
    token: Token,
    value?: TokenValue,
    className?: string,
    onSetTokenValue: (value?: TokenValue) => void,
    includeUnsetAsValue?: boolean
    highlight: boolean
}) => {

    const tokenValuesByLocalId = useSelector(selectTokenValuesByLocalId);
    const tokenValues = useSelector(selectTokenValueLocalIdsWithTokenName(token.name))
        .map(tvLocalId => tokenValuesByLocalId[tvLocalId]);

    if (includeUnsetAsValue === true) {
        tokenValues.unshift(UNSET_TOKEN_VALUE);
    }

    return (
        <SelectBase
            values={tokenValues}
            value={value}
            label={token.label}
            placeholder={token.name}
            className={className}
            onSetValue={onSetTokenValue}
            highlight={highlight}
        />
    );
};