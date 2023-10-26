import React, { useState } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import Checkbox from "@material-ui/core/Checkbox";
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
        hideAtZero: false,
    })
    counts.push({
        name: "notReady",
        description: "not ready",
        value: notReady,
        hideAtZero: true,
    })
    return counts
}


const CreateNameGroupCompose = ({ createNameGroupLocalIds, ruleName }: { createNameGroupLocalIds: string[], ruleName: string }) => {
    const [collapsed, setCollapsed] = useState<boolean>(false);

    const dispatch = useDispatch()

    const createNameGroups: CreateNameGroup[] = useSelector(selectCreateNameGroupsWithLocalIds(createNameGroupLocalIds))
    const globalState: State = useSelector(selectGlobalState)
    const selectedCreateNameGroupLocalIds: Set<string> = useSelector(selectSelectedCreateNameGroupIds)
    const renderedCompleteCreateNamesByCreateNameGroupLocalId: Record<string, string> = useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)

    const readyCounts = createReadyCounts(createNameGroups, renderedCompleteCreateNamesByCreateNameGroupLocalId)

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
            <ul className="create-name-groups-composers">
                {
                    createNameGroups.map((createNameGroup) => {
                        return (
                            <li key={createNameGroup.localId}>
                                <CreateNameGroupComposer createNameGroup={createNameGroup} includeTools />
                            </li>
                        )
                    })
                }
            </ul>
        )
    }

    return (
        <div className="create-name-groups-by-rule">
            <Grid container>
                <Grid item xs={3}>
                    <div className="create-name-groups-tools">
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
                <Grid item xs={3}>
                    <div className="create-name-groups-name">{ruleName}</div>
                </Grid>
                <Grid item xs={3}>
                    <Counts
                        counts={readyCounts}
                    />
                </Grid>
            </Grid>
            {collapsed ? <div className="composers-placeholder">...</div> : renderComposers()}
        </div>
    )
}

export default CreateNameGroupCompose;