import React, { useState, useRef } from 'react';
import { batch } from 'react-redux';
import { useAppSelector as useSelector } from '../../hooks';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import SearchIcon from '@material-ui/icons/Search';
import FormControl from '@material-ui/core/FormControl';
import DeleteOutlineOutlinedIcon from '@material-ui/icons/DeleteOutlineOutlined';
import Tooltip from '@material-ui/core/Tooltip';
import _ from 'lodash';

import { selectCreateNameGroupIdsByPrimaryRule, selectFilterStatus, selectSelectedCreateNameGroupIds } from '../../selectors/names';
import { CreateNameGroup, Rule } from '../../models';
import { RuleSelect } from '../names-create/name-composer/select';
import { selectGlobalState, selectVisibleRules } from '../../selectors/global';
import { addCreateNameGroupsToSearchCriteria, clearCreateNameGroupsFilter, clearCreateNameGroupsSelection, createNamesWithGroupForRule, deleteGroupsWithNames, setCreateNameGroupsFilterFromSearchCriteria, setCreateNameGroupsSelectionFromSearchCriteria, setSearchVisibility } from '../../actions/names';
import CreateNameGroupComposer from '../names-create/name-composer/name-composer';
import { useDispatch } from '../../utils/redux';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';



const useStyles = makeStyles((theme) => ({
    formContainer: {
        padding: '1em',
    },
    ruleSelect: {
        marginBottom: '2em',
        '&:not(:last-child)': {
            marginBottom: '1em',
        },
    },
    cngComposer: {
        margin: '0 0 2em 1em',
    },
    buttonsContainer: {
        textAlign: 'center',
        '& > div:not(:last-child)': {
            marginRight: '1em',
        },
    },
}));


interface RuleAndCreateNameGroupState {
    rule: Rule
    createNameGroup: CreateNameGroup
}


const FindAndFilterButton = ({ small, className }: { small: boolean, className?: string }) => {
    const classes = useStyles();

    const [open, setOpen] = useState<boolean>(false);
    const [ruleAndCreateNameGroup, setRuleAndCreateNameGroup] = useState<RuleAndCreateNameGroupState | undefined>();

    const dispatch = useDispatch();

    const createNameGroupIdsByPrimaryRule = useSelector(selectCreateNameGroupIdsByPrimaryRule);
    const filterEnabled = useSelector(selectFilterStatus);
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds);
    const visibleRules = useSelector(selectVisibleRules);
    const globalState = useSelector(selectGlobalState);

    const clearRuleAndCreateNameGroup = () => {
        if (ruleAndCreateNameGroup) {
            batch(() => {
                dispatch(deleteGroupsWithNames([ruleAndCreateNameGroup.createNameGroup.localId], globalState));
                setRuleAndCreateNameGroup(undefined);
            });
        }
    };

    const handleSetRule = (rule?: Rule) => {
        if (!rule) {
            if (ruleAndCreateNameGroup) {
                clearRuleAndCreateNameGroup();
                return;
            }
            setRuleAndCreateNameGroup(undefined);
            return;
        }

        const actionPayload = createNamesWithGroupForRule(rule.name, globalState, false);
        const newState: RuleAndCreateNameGroupState = {
            rule,
            createNameGroup: actionPayload.createNameGroups[0]
        };

        if (ruleAndCreateNameGroup) {
            batch(() => {
                dispatch(deleteGroupsWithNames([ruleAndCreateNameGroup.createNameGroup.localId], globalState));
                dispatch(actionPayload);
                dispatch(addCreateNameGroupsToSearchCriteria([newState.createNameGroup.localId]));
                setRuleAndCreateNameGroup(newState);
            });
        } else {
            batch(() => {
                dispatch(actionPayload);
                dispatch(addCreateNameGroupsToSearchCriteria([newState.createNameGroup.localId]));
                setRuleAndCreateNameGroup(newState);
            });
        }
    };

    const handleClickFindAndSelect = () => {
        dispatch(setCreateNameGroupsSelectionFromSearchCriteria());
    };

    const handleClickClearSelection = () => {
        batch(() => {
            dispatch(clearCreateNameGroupsSelection());
            clearRuleAndCreateNameGroup();
        });
    };

    const handleClickFilter = () => {
        dispatch(setCreateNameGroupsFilterFromSearchCriteria());
    };

    const handleClickClearFilter = () => {
        batch(() => {
            dispatch(clearCreateNameGroupsFilter());
            clearRuleAndCreateNameGroup();
        });
    };

    const renderClearSelectionButton = () => {
        const disabled = !ruleAndCreateNameGroup && selectedCreateNameGroupLocalIds.size == 0;
        const button = <Button
            onClick={handleClickClearSelection}
            color="secondary"
            size="small"
            disabled={disabled}
            aria-label='clear-selection'
        >
            <DeleteOutlineOutlinedIcon />
        </Button>;

        return disabled ? button : (
            <Tooltip title="Clear Selection">
                {button}
            </Tooltip>
        );
    };

    const renderClearFilterButton = () => {
        const disabled = !ruleAndCreateNameGroup && !filterEnabled;
        const button = <Button
            onClick={handleClickClearFilter}
            color="secondary"
            size="small"
            disabled={disabled}
            aria-label='clear-filter'
        >
            <DeleteOutlineOutlinedIcon />
        </Button>;

        return disabled ? button : (
            <Tooltip title="Clear Filter">
                {button}
            </Tooltip>
        );
    };

    const handleClick = (open: boolean) => {
        batch(() => {
            setOpen(open);
            // dispatch(setSearchVisibility(open));
        });
    };

    return (
        <ToolbarButtonWithPopper
            text="Find and Filter"
            iconComponent={<SearchIcon />}
            variant={small ? 'compact' : 'full'}
            color="#9747FF"
            popperComponent={
                <FormControl component="fieldset" id="find-and-filter-dialogue" className={classes.formContainer}>
                    <RuleSelect
                        values={visibleRules}
                        value={ruleAndCreateNameGroup?.rule}
                        label="Common Element"
                        placeholder="common element"
                        onSetRule={handleSetRule}
                        className={classes.ruleSelect}
                    />
                    {
                        ruleAndCreateNameGroup?.createNameGroup &&
                        <CreateNameGroupComposer
                            createNameGroup={ruleAndCreateNameGroup.createNameGroup}
                            includeRuleCounterIncrementer={false}
                            className={classes.cngComposer}
                            checkForDuplicates={false}
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
                                aria-label='perform-find-and-select'
                            >
                                Find + Select
                            </Button>
                            {renderClearSelectionButton()}
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
                                aria-label='perform-filter'
                            >
                                Filter
                            </Button>
                            {renderClearFilterButton()}
                        </ButtonGroup>
                    </div>
                </FormControl >
            }
            popperId="find-and-filter-dialogue"
            onClickOrPopperChange={handleClick}
            popperOpen={open}
            disabled={Object.keys(createNameGroupIdsByPrimaryRule).length == 0 && !filterEnabled}
            className={className}
        />
    );
};

export default FindAndFilterButton;