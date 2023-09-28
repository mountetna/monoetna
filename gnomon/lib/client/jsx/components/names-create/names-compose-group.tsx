import React, { useState, useRef, ReactNode } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector, useDispatch } from 'react-redux'

import { CreateName, CreateNameGroup } from "../../models";
import CreateNameGroupComposer from "./name-composer/name-composer";
import { selectCreateNameGroupByIds, selectCreateNames } from "../../selectors/names";



// move to separate module bc global toolbar uses it
const createReadyStatuses = (createNameGroups: CreateNameGroup[], createNames: Record<string, CreateName>): ReactNode => {
    const counts = { ready: 1, notReady: 1 }
    return (
        <div className="names-status">
            <div>{counts.ready}</div>
            <div>{counts.notReady}</div>
        </div>
    )
}


const CreateNameGroupCompose = ({ createNameGroupIds, ruleName }: { createNameGroupIds: string[], ruleName: string }) => {
    const createNameGroups: CreateNameGroup[] = useSelector(selectCreateNameGroupByIds(createNameGroupIds))
    const createNames: Record<string, CreateName> = useSelector(selectCreateNames)

    return (
        <div className="create-name-groups-by-rule">
            <Grid container>
                <Grid item xs={3}>
                    <div className="create-name-groups-tools">tools</div>
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