import React, { useState, useRef } from "react";
import { useSelector, useDispatch } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import AddCircleOutlineIcon from "@material-ui/icons/AddCircleOutline";

import { selectRulesByName, selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsByRuleName, selectRuleTokensByLocalId, selectTokenValueLocalIdsByTokenName } from "../../selectors/rules";
import { createNamesWithGroupForRule } from "../../actions/names";
import { Rule, RuleParent, RuleToken } from "../../models";


const useStyles = makeStyles((theme) => ({
    addNamesContainer: {
        display: "inline-block",
    },
}));


const AddNamesButton = () => {
    const classes = useStyles()
    const [open, setOpen] = useState<boolean>(false);
    const anchorEl = useRef(null)
    const dispatch = useDispatch()
    
    const rules: Record<string, Rule> = useSelector(selectRulesByName)
    const ruleParentLocalIdsByRuleName: Record<string, string[]> = useSelector(selectRuleParentLocalIdsByRuleName)
    const ruleParentsByLocalId: Record<string, RuleParent> = useSelector(selectRuleParentsByLocalId)
    const ruleTokenLocalIdsByRuleName: Record<string, string[]> = useSelector(selectRuleTokenLocalIdsByRuleName)
    const ruleTokensByLocalId: Record<string, RuleToken> = useSelector(selectRuleTokensByLocalId)
    const tokenValueLocalIdsByTokenName: Record<string, string[]> = useSelector(selectTokenValueLocalIdsByTokenName)

    const handleToggle = () => {
        setOpen(prev => !prev);
    };
    const handleClose = () => {
        setOpen(false);
    };
    const handleClickRule = (ruleName: string) => {
        dispatch(createNamesWithGroupForRule(
            ruleName,
            ruleParentLocalIdsByRuleName,
            ruleParentsByLocalId,
            ruleTokenLocalIdsByRuleName,
            ruleTokensByLocalId,
            tokenValueLocalIdsByTokenName
        ))
        handleClose();
    };

    return (
        <div className={classes.addNamesContainer}>
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
                                        Object.keys(rules).map((ruleName) =>
                                            <MenuItem
                                                onClick={() => handleClickRule(ruleName)}
                                                key={ruleName}
                                                disableRipple
                                            >
                                                {ruleName}
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