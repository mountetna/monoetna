import React, { useState, useRef } from "react";
import { useSelector } from 'react-redux'
import Button from "@material-ui/core/Button";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown';

import { CreateName, CreateNameGroup, Rule, Token, TokenValue } from "../../models";
import { selectRules, selectTokens } from "../../selectors/rules";
import { selectCreateNameById, selectCreateNamesByIds } from "../../selectors/names";




const TokenMenu = ({ token, initialValue }: { token: Token, initialValue?: TokenValue }) => {
    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null)
    const [value, setValue] = useState<TokenValue | undefined>(initialValue)

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
                startIcon={<KeyboardArrowDownIcon />}
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



const CreateNameElementsEditor = ({ createName, rule }: { createName: CreateName, rule: Rule }) => {
    const tokens: Record<string, Token> = useSelector(selectTokens)

    return (
        <span className="create-name-elements-editor">
            {
                rule.tokenNames.map((tokenName, idx) => {
                    return (
                        <React.Fragment key={idx}>
                            <TokenMenu token={tokens[tokenName]}></TokenMenu>
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