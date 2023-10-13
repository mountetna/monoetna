import React, { useState, useRef } from "react";
import { useSelector, useDispatch, batch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import ButtonGroup from "@material-ui/core/ButtonGroup";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import SearchIcon from '@material-ui/icons/Search';
import FormControl from '@material-ui/core/FormControl';
import Grid from '@material-ui/core/Grid';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import * as _ from "lodash"

import { selectCreateNameGroupIdsByPrimaryRule } from "../../selectors/names";
import { CreateNameGroup, Rule, RuleParent, RuleToken } from "../../models";
import { RuleSelect } from "../names-create/name-composer/select";
import { selectVisibleRules } from "../../selectors/global";
import { addCreateNameGroupsToSearchCriteria, clearCreateNameGroupsFilter, clearCreateNameGroupsSelection, createNamesWithGroupForRule, deleteGroupsWithNames, setCreateNameGroupsFilterFromSearchCriteria, setCreateNameGroupsSelectionFromSearchCriteria } from "../../actions/names";
import { selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsByRuleName, selectRuleTokensByLocalId, selectTokenValueLocalIdsByTokenName } from "../../selectors/rules";
import CreateNameGroupComposer from "../names-create/name-composer/name-composer";



const useStyles = makeStyles((theme) => ({
    findAndReplaceContainer: {
        display: "inline-block",
    },
}));


interface RuleAndCreateNameGroupState {
    rule: Rule
    createNameGroup: CreateNameGroup
}


const FindAndFilterButton = () => {
    const classes = useStyles()

    const [open, setOpen] = useState<boolean>(false);
    const [ruleAndCreateNameGroup, setRuleAndCreateNameGroup] = useState<RuleAndCreateNameGroupState | undefined>()

    const dispatch = useDispatch()
    const anchorEl = useRef(null)

    const createNameGroupIdsByPrimaryRule: Record<string, string[]> = useSelector(selectCreateNameGroupIdsByPrimaryRule)
    const visibleRules: Rule[] = useSelector(selectVisibleRules)
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

    const handleSetRule = (rule?: Rule) => {
        if (!rule) {
            if (ruleAndCreateNameGroup) {
                batch(() => {
                    dispatch(deleteGroupsWithNames([ruleAndCreateNameGroup.createNameGroup.localId]))
                    setRuleAndCreateNameGroup()
                })
                return
            }
            setRuleAndCreateNameGroup()
            return
        }

        const actionPayload = createNamesWithGroupForRule(
            rule.name,
            ruleParentLocalIdsByRuleName,
            ruleParentsByLocalId,
            ruleTokenLocalIdsByRuleName,
            ruleTokensByLocalId,
            tokenValueLocalIdsByTokenName,
        )
        const newState: RuleAndCreateNameGroupState = {
            rule,
            createNameGroup: actionPayload.createNameGroups[0]
        }

        if (ruleAndCreateNameGroup) {
            batch(() => {
                dispatch(deleteGroupsWithNames([ruleAndCreateNameGroup.createNameGroup.localId]))
                dispatch(actionPayload)
                dispatch(addCreateNameGroupsToSearchCriteria([newState.createNameGroup.localId]))
                setRuleAndCreateNameGroup(newState)
            })
        } else {
            batch(() => {
                dispatch(actionPayload)
                dispatch(addCreateNameGroupsToSearchCriteria([newState.createNameGroup.localId]))
                setRuleAndCreateNameGroup(newState)
            })
        }
    }

    const handleClickFindAndSelect = () => {
        dispatch(setCreateNameGroupsSelectionFromSearchCriteria())
    }

    const handleClickClearSelection = () => {
        dispatch(clearCreateNameGroupsSelection())
    }

    const handleClickFilter = () => {
        dispatch(setCreateNameGroupsFilterFromSearchCriteria())
    }

    const handleClickClearFilter = () => {
        dispatch(clearCreateNameGroupsFilter())
    }

    return (
        <div className={classes.replaceInSelectionContainer}>
            <Button
                startIcon={<SearchIcon />}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label="Find and Filter"
                aria-haspopup="true"
                aria-controls={open ? "find-and-filter-dialogue" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
                disabled={Object.keys(createNameGroupIdsByPrimaryRule).length == 0}
            >
                Find and Filter
            </Button>
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
                                <FormControl component="fieldset" id="find-and-filter-dialogue">
                                    <RuleSelect
                                        values={visibleRules}
                                        value={ruleAndCreateNameGroup?.rule}
                                        label="Rule"
                                        placeholder="rule"
                                        onSetRule={handleSetRule}
                                    />
                                    {
                                        ruleAndCreateNameGroup?.createNameGroup &&
                                        <CreateNameGroupComposer
                                            createNameGroup={ruleAndCreateNameGroup.createNameGroup}
                                        />
                                    }
                                    {/* TODO: does this need to be a grid? */}
                                    <Grid container>
                                        <Grid item xs={6}>
                                            <ButtonGroup
                                                variant="contained"
                                                color="primary"
                                                disableRipple
                                                disableElevation
                                            >
                                                <Button
                                                    onClick={handleClickFindAndSelect}
                                                    disabled={!ruleAndCreateNameGroup}
                                                >
                                                    Find + Select
                                                </Button>
                                                <Button
                                                    onClick={handleClickClearSelection}
                                                    color="secondary"
                                                    size="small"
                                                    disableElevation
                                                >
                                                    <DeleteOutlineOutlinedIcon />
                                                </Button>
                                            </ButtonGroup>
                                        </Grid>
                                        <Grid item xs={6}>
                                            <ButtonGroup
                                                variant="contained"
                                                color="primary"
                                                disableRipple
                                                disableElevation
                                            >
                                                <Button
                                                    onClick={handleClickFilter}
                                                    disabled={!ruleAndCreateNameGroup}
                                                >
                                                    Filter
                                                </Button>
                                                <Button
                                                    onClick={handleClickClearFilter}
                                                    color="secondary"
                                                    size="small"
                                                    disableElevation
                                                >
                                                    <DeleteOutlineOutlinedIcon />
                                                </Button>
                                            </ButtonGroup>
                                        </Grid>
                                    </Grid>
                                </FormControl>
                            </ClickAwayListener>
                        </Paper>
                    </Grow>
                )}
            </Popper>
        </div>
    )
}

export default FindAndFilterButton