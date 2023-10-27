import React, { useState } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import Checkbox from "@material-ui/core/Checkbox";
import { makeStyles } from '@material-ui/core/styles';
import * as _ from "lodash"
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import AddCircleOutlineIcon from "@material-ui/icons/AddCircleOutline";

import { CreateNameGroup } from "../../models";
import CreateNameGroupComposer from "./name-composer/name-composer";
import { selectCreateNameGroupsWithLocalIds, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId, selectSelectedCreateNameGroupIds } from "../../selectors/names";
import { deleteGroupsWithNames, createNamesWithGroupForRule, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection } from "../../actions/names";
import { useDispatch } from "../../utils/redux";
import { State } from "../../store";
import { selectGlobalState } from "../../selectors/global";
import Counts, { Count } from "./counts";



function createReadyCounts(
    createNameGroups: CreateNameGroup[],
    renderedCompleteCreateNamesByCreateNameGroupLocalId: Record<string, string>,
): Count[] {

    let total = 0
    let notReady = 0

    for (const cng of createNameGroups) {
        total += 1
        if (!(cng.localId in renderedCompleteCreateNamesByCreateNameGroupLocalId)) {
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
    container: {
        display: "inline-block",
        minWidth: "30em",
    },
    headerContainer: {
        alignItems: "center",
    },
    tools: {
        textAlign: "left",
    },
    ruleNameTitle: {
        textAlign: "center",
    },
    readyCounts: {
        textAlign: "right",
        "& .count-not-ready": {
            color: "red",
        },
    },
    composersContainer: {
        maxHeight: "80vh",
        overflow: "scroll",
        border: "1px solid black",
        borderRadius: "6px",
        "&.some-not-ready": {
            borderColor: "red"
        },
        padding: "2em 1em",
        "&.collapsed": {
            transition: "background 0.2s ease-out",
        },
        "&.collapsed:hover": {
            cursor: "pointer",
            background: "rgba(0, 0, 0, 0.05)",
        },
    },
    placeholder: {

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

    const createNameGroups: CreateNameGroup[] = useSelector(selectCreateNameGroupsWithLocalIds(createNameGroupLocalIds))
    const globalState: State = useSelector(selectGlobalState)
    const selectedCreateNameGroupLocalIds: Set<string> = useSelector(selectSelectedCreateNameGroupIds)
    const renderedCompleteCreateNamesByCreateNameGroupLocalId: Record<string, string> = useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)

    const readyCounts = createReadyCounts(createNameGroups, renderedCompleteCreateNamesByCreateNameGroupLocalId)
    const someNotReady = readyCounts[1].value > 0

    const handleClickAdd = () => {
        dispatch(createNamesWithGroupForRule(ruleName, globalState, true))
    }

    const handleClickSelect = (event: React.ChangeEvent) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection(createNameGroupLocalIds))
            return
        }
        dispatch(removeCreateNameGroupsFromSelection(createNameGroupLocalIds))
    }

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames(createNameGroupLocalIds, globalState))
    }

    const renderComposers = (): React.ReactNode => {
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
        <div className={classes.container}>
            <Grid container className={classes.headerContainer}>
                <Grid item xs={4}>
                    <div className={classes.tools}>
                        <Checkbox
                            checked={_.every(createNameGroups, (cng: CreateNameGroup) => selectedCreateNameGroupLocalIds.has(cng.localId))}
                            onChange={handleClickSelect}
                            inputProps={{ 'aria-label': 'Select the Name Group' }}
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
                {collapsed
                    ? (
                        <div className={classes.placeholder}>...</div>
                    )
                    : renderComposers()}
            </div>
        </div>
    )
}

export default CreateNameGroupCompose;