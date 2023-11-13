import React, { useState } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import Checkbox from "@material-ui/core/Checkbox";
import { makeStyles } from '@material-ui/core/styles';
import _ from "lodash";
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import AddCircleOutlineIcon from "@material-ui/icons/AddCircleOutline";
import { useVirtualizer } from '@tanstack/react-virtual'

import { CreateNameGroup } from "../../models";
import CreateNameGroupComposer from "./name-composer/name-composer";
import { selectComposeErrorsByCreateNameGroupLocalId, selectCreateNameGroupsWithLocalIds, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId, selectSelectedCreateNameGroupIds } from "../../selectors/names";
import { deleteGroupsWithNames, createNamesWithGroupForRule, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection } from "../../actions/names";
import { useDispatch } from "../../utils/redux";
import { selectGlobalState } from "../../selectors/global";
import Counts, { Count } from "./counts";



function createReadyCounts(
    createNameGroups: CreateNameGroup[],
    renderedCompleteCreateNamesByCreateNameGroupLocalId: Record<string, string>,
    composeErrorsByCreateNameGroupLocalId: Record<string, boolean>,
): Count[] {

    let total = 0
    let notReady = 0

    for (const cng of createNameGroups) {
        total += 1
        if (
            !(cng.localId in renderedCompleteCreateNamesByCreateNameGroupLocalId)
            || composeErrorsByCreateNameGroupLocalId[cng.localId]
        ) {
            notReady += 1
        }
    }

    const counts: Count[] = []
    counts.push({
        name: "total",
        description: `name${total > 1 ? "s" : ""}`,
        value: total,
    })
    counts.push({
        name: "not-ready",
        description: "not ready",
        value: notReady,
        hideAtZero: true,
    })
    return counts
}


const useStyles = makeStyles((theme) => ({
    gridContainer: {
        display: "inline-block",
        width: "34em",
        maxWidth: "100%",
        marginBottom: "4em",
    },
    innerContainer: {},
    headerContainer: {
        alignItems: "center",
        marginBottom: "0.5em",
    },
    tools: {
        textAlign: "left",
    },
    checkbox: {
        padding: "0",
    },
    ruleNameTitle: {
        textAlign: "center",
        fontWeight: "bold",
        fontSize: "18px",
    },
    readyCounts: {
        fontWeight: "bold",
        textAlign: "right",
        "& .count-not-ready": {
            color: "red",
        },
    },
    composersContainer: {
        maxHeight: "50vh",
        overflow: "scroll",
        border: "2px solid black",
        borderRadius: "6px",
        "&.some-not-ready": {
            borderColor: "red"
        },
        padding: "2em",
        "& .placeholder": {
            display: "none",
        },
        "&.collapsed": {
            transition: "background 0.2s ease-out",
            "& > .placeholder": {
                display: "block",
            },
            "& > :not(.placeholder)": {
                display: "none",
            },
        },
        "&.collapsed:hover": {
            cursor: "pointer",
            background: "rgba(0, 0, 0, 0.05)",
        },
    },
    composerList: {

    },
    composer: {
        textAlign: "left",
    },
}));



interface ListItemData {
    className: string
    createNameGroups: CreateNameGroup[]
}


const ListItem = ({ data, index, style }: { data: ListItemData, index: number, style: object }) => {
    return (
        <div
            style={style}
            className={data.className}
        >
            <CreateNameGroupComposer createNameGroup={data.createNameGroups[index]} includeTools />
        </div>
    )
}


const CreateNameGroupRuleGroupCompose = ({ createNameGroupLocalIds, ruleName }: { createNameGroupLocalIds: string[], ruleName: string }) => {
    const [collapsed, setCollapsed] = useState<boolean>(false);

    const dispatch = useDispatch()
    const classes = useStyles()

    const createNameGroups = useSelector(selectCreateNameGroupsWithLocalIds(createNameGroupLocalIds))
    const globalState = useSelector(selectGlobalState)
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds)
    const renderedCompleteCreateNamesByCreateNameGroupLocalId = useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)
    const composeErrorsByCreateNameGroupLocalId = useSelector(selectComposeErrorsByCreateNameGroupLocalId)

    const readyCounts = createReadyCounts(
        createNameGroups,
        renderedCompleteCreateNamesByCreateNameGroupLocalId,
        composeErrorsByCreateNameGroupLocalId
    )
    const someNotReady = readyCounts[1].value > 0

    // handle rendering composers
    const virtualizerParentRef = React.useRef<HTMLDivElement>(null)
    const virtualizer = useVirtualizer({
        count: createNameGroups.length,
        getScrollElement: () => virtualizerParentRef.current,
        estimateSize: () => 65,
        overscan: 20,
    })

    const renderComposers = () => {
        return (
            <div
                className={classes.composerList}
                style={{
                    height: virtualizer.getTotalSize(),
                    width: "100%",
                    position: "relative",
                }}
            >
                {
                    virtualizer.getVirtualItems().map((row) => {
                        const createNameGroup = createNameGroups[row.index]

                        return (
                            <div
                                key={createNameGroup.localId}
                                data-index={row.index}
                                className={classes.composer}
                                style={{
                                    position: "absolute",
                                    top: 0,
                                    left: 0,
                                    width: "100%",
                                    height: `${row.size}px`,
                                    transform: `translateY(${row.start}px)`,
                                }}
                            >
                                <div ref={virtualizer.measureElement}>
                                    <CreateNameGroupComposer createNameGroup={createNameGroup} includeTools />
                                </div>
                            </div>
                        )
                    })
                }
            </div>
        )
    }

    const handleClickAdd = () => {
        dispatch(createNamesWithGroupForRule(ruleName, globalState, true))
    }

    const handleClickSelect = (event: React.ChangeEvent<HTMLInputElement>) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection(createNameGroupLocalIds))
            return
        }
        dispatch(removeCreateNameGroupsFromSelection(createNameGroupLocalIds))
    }

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames(createNameGroupLocalIds, globalState))
    }

    return (
        <div className={classes.gridContainer}>
            <div className={classes.innerContainer}>
                <Grid container className={classes.headerContainer}>
                    <Grid item xs={4}>
                        <div className={classes.tools}>
                            <Checkbox
                                checked={_.every(createNameGroups, (cng) => selectedCreateNameGroupLocalIds.has(cng.localId))}
                                onChange={handleClickSelect}
                                inputProps={{ 'aria-label': 'Select the Name Group' }}
                                className={classes.checkbox}
                            />
                            <ButtonBase
                                onClick={() => setCollapsed(prev => !prev)}
                                aria-label={"Toggle Expand/Collapse"}
                                disableRipple
                                disableTouchRipple
                            >
                                {collapsed ? <UnfoldMoreOutlinedIcon /> : <UnfoldLessOutlinedIcon />}
                            </ButtonBase>
                            <ButtonBase
                                onClick={handleClickAdd}
                                aria-label="Add Name"
                                disableRipple
                                disableTouchRipple
                            >
                                <AddCircleOutlineIcon />
                            </ButtonBase>
                            <ButtonBase
                                onClick={handleClickDelete}
                                aria-label="Delete Name"
                                disableRipple
                                disableTouchRipple
                            >
                                <DeleteOutlineOutlinedIcon />
                            </ButtonBase>
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
                    className={`${classes.composersContainer} ${someNotReady ? "some-not-ready" : ""} ${collapsed ? "collapsed" : ""}`}
                    onClick={() => collapsed && setCollapsed(false)}
                    ref={virtualizerParentRef}
                >
                    <div className="placeholder">...</div>
                    {renderComposers()}
                </div>
            </div>
        </div>
    )
}

export default CreateNameGroupRuleGroupCompose;