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

import { CreateName, CreateNameGroup, Rule, TOKEN_VALUE_PLACEHOLDER, Token, TokenValue } from "../../models";
import { selectRules, selectTokens } from "../../selectors/rules";
import { selectCreateNameById, selectCreateNamesByIds } from "../../selectors/names";
import { setCreateNameTokenValue } from "../../actions/names";


const useStyles = makeStyles((theme) => ({
    selectTokenValueContainer: {
        display: "inline-block",
    },
    currentTokenValue: {
        fontWeight: "bold"
    },
    currentTokenValueUnset: {
        color: "red"
    }
}));


const TokenMenu = ({ token, value, handleSetTokenValue }:
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
                    className={classes.currentTokenValue + " " + (!value ? `${classes.currentTokenValueUnset}` : "")}
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
                                    <div
                                        key={token.label}
                                    >
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



const CreateNameElementsEditor = ({ createName, rule }: { createName: CreateName, rule: Rule }) => {
    const tokens: Record<string, Token> = useSelector(selectTokens)
    const dispatch = useDispatch()

    const setTokenValue = (tokenValue: TokenValue, tokenIdx: number) => {
        dispatch(setCreateNameTokenValue(createName.localId, tokenValue, tokenIdx))
    }

    return (
        <span className="create-name-elements-editor">
            {
                rule.tokenNames.map((tokenName, idx) => {
                    const currentValue = createName.tokenValues[idx]

                    return (
                        <React.Fragment key={idx}>
                            <TokenMenu
                                token={tokens[tokenName]}
                                value={currentValue != TOKEN_VALUE_PLACEHOLDER ? currentValue : undefined}
                                handleSetTokenValue={(tokenValue: TokenValue) => setTokenValue(tokenValue, idx)}
                            ></TokenMenu>
                        </React.Fragment>
                    )
                })
            }
        </span>
    )
}



const CreateNameGroupComposer = ({ createNameGroup }: { createNameGroup: CreateNameGroup }) => {
    const createNames: CreateName[] = useSelector(selectCreateNamesByIds(createNameGroup.createNameIds))
    const primaryCreateName: CreateName = useSelector(selectCreateNameById(createNameGroup.primaryCreateNameId))
    const allRules: Record<string, Rule> = useSelector(selectRules)
    const orderedRuleNames: string[] = []

    const ruleNamesToScan: string[] = [primaryCreateName.ruleName]
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop()
        if (!ruleName) {
            break
        }
        const rule = allRules[ruleName]

        orderedRuleNames.unshift(rule.name)

        if (rule.parentRuleNames.length) {
            ruleNamesToScan.push(...rule.parentRuleNames)
        }
    }

    return (
        <div className="create-name-group-composer">
            <span className="create-name-group-composer-tools"></span>
            {
                orderedRuleNames.map((ruleName) => {
                    const createName = createNames.find((createName) => createName.ruleName == ruleName)
                    if (!createName) {
                        console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`)
                        return
                    }
                    const rule = allRules[ruleName]

                    return (
                        <React.Fragment key={rule.name}>
                            <CreateNameElementsEditor createName={createName} rule={rule} />
                        </React.Fragment>
                    )
                })
            }
        </div>
    )
}


export default CreateNameGroupComposer;