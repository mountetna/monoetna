import React, { useState, useRef } from "react";
import { useSelector, batch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import ButtonGroup from "@material-ui/core/ButtonGroup";
import SearchIcon from '@material-ui/icons/Search';
import FormControl from '@material-ui/core/FormControl';
import Grid from '@material-ui/core/Grid';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import _ from "lodash"

import { selectCreateNameGroupIdsByPrimaryRule, selectFilterStatus, selectSelectedCreateNameGroupIds } from "../../selectors/names";
import { CreateNameGroup, Rule } from "../../models";
import { RuleSelect } from "../names-create/name-composer/select";
import { selectGlobalState, selectVisibleRules } from "../../selectors/global";
import { addCreateNameGroupsToSearchCriteria, clearCreateNameGroupsFilter, clearCreateNameGroupsSelection, createNamesWithGroupForRule, deleteGroupsWithNames, setCreateNameGroupsFilterFromSearchCriteria, setCreateNameGroupsSelectionFromSearchCriteria } from "../../actions/names";
import CreateNameGroupComposer from "../names-create/name-composer/name-composer";
import { useDispatch } from "../../utils/redux";
import ToolbarButtonWithPopper from "./toolbar-button-with-popper";



const useStyles = makeStyles((theme) => ({
    container: {
    },
    buttonsContainer: {

    },
}));


interface RuleAndCreateNameGroupState {
    rule: Rule
    createNameGroup: CreateNameGroup
}


const FindAndFilterButton = ({ small }: { small: boolean }) => {
    const classes = useStyles()

    const [open, setOpen] = useState<boolean>(false);
    const [ruleAndCreateNameGroup, setRuleAndCreateNameGroup] = useState<RuleAndCreateNameGroupState | undefined>()

    const dispatch = useDispatch()

    const createNameGroupIdsByPrimaryRule = useSelector(selectCreateNameGroupIdsByPrimaryRule)
    const filterEnabled = useSelector(selectFilterStatus)
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds)
    const visibleRules = useSelector(selectVisibleRules)
    const globalState = useSelector(selectGlobalState)

    const handleSetRule = (rule?: Rule) => {
        if (!rule) {
            if (ruleAndCreateNameGroup) {
                batch(() => {
                    dispatch(deleteGroupsWithNames([ruleAndCreateNameGroup.createNameGroup.localId], globalState))
                    setRuleAndCreateNameGroup(undefined)
                })
                return
            }
            setRuleAndCreateNameGroup(undefined)
            return
        }

        const actionPayload = createNamesWithGroupForRule(rule.name, globalState, false)
        const newState: RuleAndCreateNameGroupState = {
            rule,
            createNameGroup: actionPayload.createNameGroups[0]
        }

        if (ruleAndCreateNameGroup) {
            batch(() => {
                dispatch(deleteGroupsWithNames([ruleAndCreateNameGroup.createNameGroup.localId], globalState))
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
        batch(() => {
            dispatch(clearCreateNameGroupsSelection())
            setRuleAndCreateNameGroup(undefined)
        })
    }

    const handleClickFilter = () => {
        dispatch(setCreateNameGroupsFilterFromSearchCriteria())
    }

    const handleClickClearFilter = () => {
        batch(() => {
            dispatch(clearCreateNameGroupsFilter())
            setRuleAndCreateNameGroup(undefined)
        })
    }

    return (
        <ToolbarButtonWithPopper
            text="Find and Filter"
            iconComponent={<SearchIcon />}
            variant={small ? "compact" : "full"}
            color="primary"
            popperComponent={
                <FormControl component="fieldset" id="find-and-filter-dialogue" className={classes.container}>
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
                    <div className={classes.buttonsContainer}>
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
                                disabled={!ruleAndCreateNameGroup && selectedCreateNameGroupLocalIds.size == 0}
                            >
                                <DeleteOutlineOutlinedIcon />
                            </Button>
                        </ButtonGroup>
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
                                disabled={!ruleAndCreateNameGroup && !filterEnabled}
                            >
                                <DeleteOutlineOutlinedIcon />
                            </Button>
                        </ButtonGroup>
                    </div>
                </FormControl >
            }
            popperId="find-and-filter-dialogue"
            onClickOrPopperChange={(open: boolean) => setOpen(open)}
            popperOpen={open}
            disabled={Object.keys(createNameGroupIdsByPrimaryRule).length == 0 && !filterEnabled}
        />
    )
}

export default FindAndFilterButton