import React from 'react';
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown';
import ClickAwayListener from '@material-ui/core/ClickAwayListener';
import Grow from '@material-ui/core/Grow';
import Paper from '@material-ui/core/Paper';
import Popper from '@material-ui/core/Popper';
import MenuItem from '@material-ui/core/MenuItem';
import MenuList from '@material-ui/core/MenuList';
import _ from 'lodash';

import { createLocalId } from '../utils/models';



const MultiOptionButton = <Option extends string>({ options, optionPrefix, onClick, ...buttonProps }: {
    options: Option[],
    optionPrefix?: string,
    onClick: (option: Option) => void,
    [buttonProps: string]: any,
}) => {
    const formId = createLocalId();

    const [open, setOpen] = React.useState(false);
    const anchorRef = React.useRef<HTMLButtonElement>(null);
    const [selectedIndex, setSelectedIndex] = React.useState(0);

    const handleToggle = () => {
        setOpen((prevOpen) => !prevOpen);
    };

    const handleClose = (event: React.MouseEvent<Document, MouseEvent>) => {
        if (anchorRef.current && anchorRef.current.contains(event.target as HTMLElement)) {
            return;
        }
        setOpen(false);
    };

    const handleClick = () => {
        onClick(options[selectedIndex]);
    };

    const handleMenuItemClick = (
        event: React.MouseEvent<HTMLLIElement, MouseEvent>,
        index: number,
    ) => {
        setSelectedIndex(index);
        setOpen(false);
    };

    return (
        <React.Fragment>
            <ButtonGroup
                variant="contained"
                aria-label="split button"
                {...buttonProps}
            >
                <Button
                    onClick={handleClick}
                    ref={anchorRef}
                    {..._.omit(buttonProps, ['className'])}
                >
                    {optionPrefix
                        ? `${optionPrefix} ${options[selectedIndex]}`
                        : options[selectedIndex]}
                </Button>
                <Button
                    size="small"
                    aria-controls={open ? formId : undefined}
                    aria-expanded={open ? 'true' : undefined}
                    aria-label="select option"
                    aria-haspopup="menu"
                    onClick={handleToggle}
                    className={`menuToggle ${buttonProps.className ? buttonProps.className : ''}`}
                    {..._.omit(buttonProps, ['className'])}
                >
                    <ArrowDropDownIcon />
                </Button>
            </ButtonGroup>
            <Popper
                open={open}
                anchorEl={anchorRef.current}
                role={undefined}
                transition
                disablePortal
                placement='bottom'
            >
                {({ TransitionProps }) => (
                    <Grow
                        {...TransitionProps}
                        style={{ transformOrigin: 'center top' }}
                    >
                        <Paper>
                            <ClickAwayListener onClickAway={handleClose}>
                                <MenuList id={formId}>
                                    {options.map((option, index) => (
                                        <MenuItem
                                            key={option}
                                            selected={index === selectedIndex}
                                            onClick={(event) => handleMenuItemClick(event, index)}
                                            disableRipple
                                        >
                                            {option}
                                        </MenuItem>
                                    ))}
                                </MenuList>
                            </ClickAwayListener>
                        </Paper>
                    </Grow>
                )}
            </Popper>
        </React.Fragment>
    );
};


export default MultiOptionButton;