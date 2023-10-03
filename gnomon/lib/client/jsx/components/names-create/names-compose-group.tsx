import React, { useState, useRef, ReactNode } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector, useDispatch } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import Checkbox from "@material-ui/core/Checkbox";
import * as _ from "lodash"
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import AddCircleOutlineIcon from "@material-ui/icons/AddCircleOutline";

import { CreateName, CreateNameGroup, RuleParent, RuleToken } from "../../models";
import CreateNameGroupComposer from "./name-composer/name-composer";
import { selectCreateNameGroupsWithLocalIds, selectCreateNamesByLocalId } from "../../selectors/names";
import { setCreateNameGroupsSelected, deleteGroupsWithNames, createNamesWithGroupForRule } from "../../actions/names";
import { selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsByRuleName, selectRuleTokensByLocalId, selectTokenValuesNamesByTokenName } from "../../selectors/rules";



// move to separate module bc global toolbar uses it
const createReadyStatuses = (createNameGroups: CreateNameGroup[], createNames: Record<string, CreateName>): ReactNode => {
    const counts = { ready: 1, notReady: 1 }
    return (
        <div className="names-status">
            <div>{counts.ready} ready</div>
            {counts.notReady && <div>{counts.notReady} not ready</div>}
        </div>
    )
}


const CreateNameGroupCompose = ({ createNameGroupIds, ruleName }: { createNameGroupIds: string[], ruleName: string }) => {
    const dispatch = useDispatch()
    const createNameGroups: CreateNameGroup[] = useSelector(selectCreateNameGroupsWithLocalIds(createNameGroupIds))
    const createNames: Record<string, CreateName> = useSelector(selectCreateNamesByLocalId)
    const ruleParentLocalIdsByRuleName: Record<string, string[]> = useSelector(selectRuleParentLocalIdsByRuleName)
    const ruleParentsByLocalId: Record<string, RuleParent> = useSelector(selectRuleParentsByLocalId)
    const ruleTokenLocalIdsByRuleName: Record<string, string[]> = useSelector(selectRuleTokenLocalIdsByRuleName)
    const ruleTokensByLocalId: Record<string, RuleToken> = useSelector(selectRuleTokensByLocalId)
    const tokenValueNamesByTokenName: Record<string, string[]> = useSelector(selectTokenValuesNamesByTokenName)

    const [collapsed, setCollapsed] = useState<Boolean>(false);

    const handleClickAdd = () => {
        dispatch(createNamesWithGroupForRule(
            ruleName,
            ruleParentLocalIdsByRuleName,
            ruleParentsByLocalId,
            ruleTokenLocalIdsByRuleName,
            ruleTokensByLocalId,
            tokenValueNamesByTokenName
        ))
    }

    const handleClickSelect = (event: React.ChangeEvent) => {
        dispatch(setCreateNameGroupsSelected(createNameGroupIds, event.target.checked))
    }

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames(createNameGroupIds))
    }

    return (
        <div className="create-name-groups-by-rule">
            <Grid container>
                <Grid item xs={3}>
                    <div className="create-name-groups-tools">
                        <Checkbox
                            checked={_.every(createNameGroups, (cng) => cng.selected)}
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
                    {createReadyStatuses(createNameGroups, createNames)}
                </Grid>
            </Grid>
            <ul className="create-name-groups-composers">
                {
                    createNameGroups.map((createNameGroup) => {
                        return (
                            <li key={createNameGroup.localId}>
                                <CreateNameGroupComposer createNameGroup={createNameGroup} />
                            </li>
                        )
                    })
                }
            </ul>
        </div>
    )
}

export default CreateNameGroupCompose;