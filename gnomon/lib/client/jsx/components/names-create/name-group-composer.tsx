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
        marginBottom: "1em",
        "&:last-child": {
            marginBottom: "0",
        },
    },
}));


const CreateNameGroupCompose = ({ createNameGroupLocalIds, ruleName }: { createNameGroupLocalIds: string[], ruleName: string }) => {
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

    const renderComposers = () => {
        return (
            <div className={classes.composerList}>
                {
                    createNameGroups.map((createNameGroup) => {
                        return (
                            <div key={createNameGroup.localId} className={classes.composer}>
                                <CreateNameGroupComposer createNameGroup={createNameGroup} includeTools />
                            </div>
                        )
                    })
                }
            </div>
        )
    }

    return (
        <div className={classes.gridContainer}>
            <div className={classes.innerContainer}>
                <Grid container className={classes.headerContainer}>
                    <Grid item xs={4}>
                        <div className={classes.tools}>
                            <Checkbox
                                checked={_.every(createNameGroups, (cng: CreateNameGroup) => selectedCreateNameGroupLocalIds.has(cng.localId))}
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
                >
                    <div className="placeholder">...</div>
                    {renderComposers()}
                </div>
            </div>
        </div>
    )
}

export default CreateNameGroupCompose;