import React, { useState, useRef } from 'react';
import Grid from '@material-ui/core/Grid';
import { batch } from 'react-redux';
import { useAppSelector as useSelector } from '../../hooks';
import ButtonBase from '@material-ui/core/ButtonBase';
import DeleteOutlineOutlinedIcon from '@material-ui/icons/DeleteOutlineOutlined';
import Checkbox from '@material-ui/core/Checkbox';
import { makeStyles } from '@material-ui/core/styles';
import _ from 'lodash';
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import AddCircleOutlineIcon from '@material-ui/icons/AddCircleOutline';
import Tooltip from '@material-ui/core/Tooltip';
import { useVirtualizer } from '@tanstack/react-virtual';

import { CreateNameGroup } from '../../models';
import CreateNameGroupComposer from './name-composer/name-composer';
import { selectComposeErrorsByCreateNameGroupLocalId, selectCreateNameGroupsWithLocalIds, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId, selectSelectedCreateNameGroupIds } from '../../selectors/names';
import { deleteGroupsWithNames, createNamesWithGroupForRule, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection } from '../../actions/names';
import { useDispatch } from '../../utils/redux';
import { selectGlobalState } from '../../selectors/global';
import Counts, { Count } from './counts';
import ConfirmationPopper from '../confirmation-popper';



function createReadyCounts(
    createNameGroups: CreateNameGroup[],
    renderedCompleteCreateNamesByCreateNameGroupLocalId: Record<string, string>,
    composeErrorsByCreateNameGroupLocalId: Record<string, boolean>,
): Count[] {

    let total = 0;
    let notReady = 0;

    for (const cng of createNameGroups) {
        total += 1;
        if (
            !(cng.localId in renderedCompleteCreateNamesByCreateNameGroupLocalId)
            || !(cng.localId in composeErrorsByCreateNameGroupLocalId)
            || composeErrorsByCreateNameGroupLocalId[cng.localId]
        ) {
            notReady += 1;
        }
    }

    const counts: Count[] = [];
    counts.push({
        name: 'total',
        description: `name${total > 1 ? 's' : ''}`,
        value: total,
    });
    counts.push({
        name: 'not-ready',
        description: 'not ready',
        value: notReady,
        hideAtZero: true,
    });
    return counts;
}


const useStyles = makeStyles((theme) => ({
    gridContainer: {
        display: 'inline-block',
        width: '34em',
        maxWidth: '100%',
        marginBottom: '4em',
    },
    innerContainer: {},
    headerContainer: {
        alignItems: 'center',
        marginBottom: '0.5em',
    },
    tools: {
        textAlign: 'left',
    },
    checkbox: {
        padding: '0',
        '&.MuiCheckbox-indeterminate': {
            color: theme.palette.secondary.main,
        },
    },
    ruleNameTitle: {
        textAlign: 'center',
        fontWeight: 'bold',
        fontSize: '18px',
    },
    readyCounts: {
        fontWeight: 'bold',
        textAlign: 'right',
        '& .count-not-ready': {
            color: 'red',
        },
    },
    composersContainer: {
        maxHeight: '50vh',
        overflow: 'auto',
        border: '2px solid black',
        borderRadius: '6px',
        '&.some-not-ready': {
            borderColor: 'red'
        },
        padding: '2em',
        '& .placeholder': {
            display: 'none',
        },
        '&.collapsed': {
            transition: 'background 0.2s ease-out',
            '& > .placeholder': {
                display: 'block',
            },
            '& > :not(.placeholder)': {
                display: 'none',
            },
        },
        '&.collapsed:hover': {
            cursor: 'pointer',
            background: 'rgba(0, 0, 0, 0.05)',
        },
    },
    composerList: {},
    composerRowContainer: {
        '&:not(:last-child) > div': {
            paddingBottom: '1em',
        }
    },
    composerContainer: {
        textAlign: 'left',
    },
}));


const CreateNameGroupRuleGroupCompose = ({ createNameGroupLocalIds, ruleName }: { createNameGroupLocalIds: string[], ruleName: string }) => {
    const [collapsed, setCollapsed] = useState<boolean>(false);
    const [deleteConfirmationPopperOpen, setDeleteConfirmationPopperOpen] = useState<boolean>(false);

    const deleteButtonRef = useRef(null);

    const dispatch = useDispatch();
    const classes = useStyles();

    const createNameGroups = useSelector(selectCreateNameGroupsWithLocalIds(createNameGroupLocalIds));
    const globalState = useSelector(selectGlobalState);
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds);
    const renderedCompleteCreateNamesByCreateNameGroupLocalId = useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId);
    const composeErrorsByCreateNameGroupLocalId = useSelector(selectComposeErrorsByCreateNameGroupLocalId);

    const readyCounts = createReadyCounts(
        createNameGroups,
        renderedCompleteCreateNamesByCreateNameGroupLocalId,
        composeErrorsByCreateNameGroupLocalId
    );
    const someNotReady = readyCounts[1].value > 0;

    // handle rendering composers with virtual
    const virtualizerParentRef = React.useRef<HTMLDivElement>(null);
    const virtualizer = useVirtualizer({
        count: createNameGroups.length,
        getScrollElement: () => virtualizerParentRef.current,
        estimateSize: () => 65,
        overscan: 500,
    });

    // TODO: switch to this rendering strategy when it's performant
    const renderComposersWithVirtual = () => {
        return (
            <div
                className={classes.composerList}
                style={{
                    height: virtualizer.getTotalSize(),
                    width: '100%',
                    position: 'relative',
                }}
            >
                {
                    virtualizer.getVirtualItems().map((row) => {
                        const createNameGroup = createNameGroups[row.index];

                        return (
                            <div
                                key={createNameGroup.localId}
                                style={{
                                    position: 'absolute',
                                    top: 0,
                                    left: 0,
                                    width: '100%',
                                    height: `${row.size}px`,
                                    transform: `translateY(${row.start}px)`,
                                }}
                                className={classes.composerRowContainer}
                            >
                                <div
                                    ref={virtualizer.measureElement}
                                    data-index={row.index}
                                    className={classes.composerContainer}
                                >
                                    <CreateNameGroupComposer createNameGroup={createNameGroup} includeTools />
                                </div>
                            </div>
                        );
                    })
                }
            </div>
        );
    };

    const renderComposers = () => {
        return (
            <div className={classes.composerList}>
                {
                    createNameGroups.map((cng) => {

                        return (
                            <div
                                key={cng.localId}
                                className={classes.composerRowContainer}
                            >
                                <div className={classes.composerContainer}>
                                    <CreateNameGroupComposer createNameGroup={cng} includeTools />
                                </div>
                            </div>
                        );
                    })
                }
            </div>
        );
    };

    const handleClickAdd = () => {
        dispatch(createNamesWithGroupForRule(ruleName, globalState, true));
    };

    const handleClickSelect = (event: React.ChangeEvent<HTMLInputElement>) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection(createNameGroupLocalIds));
            return;
        }
        dispatch(removeCreateNameGroupsFromSelection(createNameGroupLocalIds));
    };

    const handleConfirmDelete = (confirmed: boolean) => {
        if (!confirmed) {
            setDeleteConfirmationPopperOpen(false);
            return;
        }

        batch(() => {
            dispatch(deleteGroupsWithNames(createNameGroupLocalIds, globalState));
            setDeleteConfirmationPopperOpen(false);
        });
    };

    const allSelected = _.every(createNameGroups, (cng) => selectedCreateNameGroupLocalIds.has(cng.localId));
    const someSelected = _.some(createNameGroups, (cng) => selectedCreateNameGroupLocalIds.has(cng.localId));

    return (
        <div className={classes.gridContainer}>
            <div className={classes.innerContainer}>
                <Grid container className={classes.headerContainer}>
                    <Grid item xs={4}>
                        <div className={classes.tools}>
                            <Checkbox
                                checked={allSelected}
                                indeterminate={allSelected ? false : someSelected}
                                onChange={handleClickSelect}
                                inputProps={{ 'aria-label': 'Select the Name Group' }}
                                className={classes.checkbox}
                            />
                            <Tooltip title={collapsed ? 'Expand' : 'Collapse'} placement='top'>
                                <ButtonBase
                                    onClick={() => setCollapsed(prev => !prev)}
                                    aria-label={'Toggle Expand/Collapse'}
                                    disableRipple
                                    disableTouchRipple
                                >
                                    {collapsed ? <UnfoldMoreOutlinedIcon /> : <UnfoldLessOutlinedIcon />}
                                </ButtonBase>
                            </Tooltip>
                            <Tooltip title='Add Name' placement='top'>
                                <ButtonBase
                                    onClick={handleClickAdd}
                                    aria-label="Add Name"
                                    disableRipple
                                    disableTouchRipple
                                >
                                    <AddCircleOutlineIcon />
                                </ButtonBase>
                            </Tooltip>
                            <Tooltip title='Delete Names' placement='top'>
                                <ButtonBase
                                    onClick={() => setDeleteConfirmationPopperOpen(open => !open)}
                                    aria-label="Delete Names"
                                    disableRipple
                                    disableTouchRipple
                                    ref={deleteButtonRef}
                                >
                                    <DeleteOutlineOutlinedIcon />
                                </ButtonBase>
                            </Tooltip>
                            <ConfirmationPopper
                                open={deleteConfirmationPopperOpen}
                                onConfirm={handleConfirmDelete}
                                onClose={() => setDeleteConfirmationPopperOpen(false)}
                                anchorRef={deleteButtonRef}
                            />
                        </div>
                    </Grid>
                    <Grid item xs={4}>
                        <div className={classes.ruleNameTitle}>{ruleName}</div>
                    </Grid>
                    <Grid item xs={4}>
                        <Counts
                            counts={readyCounts}
                            className={classes.readyCounts}
                        />
                    </Grid>
                </Grid>
                <div
                    className={`${classes.composersContainer} ${someNotReady ? 'some-not-ready' : ''} ${collapsed ? 'collapsed' : ''}`}
                    onClick={() => collapsed && setCollapsed(false)}
                    ref={virtualizerParentRef}
                >
                    <div className="placeholder">...</div>
                    {renderComposers()}
                </div>
            </div>
        </div>
    );
};

export default CreateNameGroupRuleGroupCompose;