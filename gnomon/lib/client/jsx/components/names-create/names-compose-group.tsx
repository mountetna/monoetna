import React, { useState, useRef, ReactNode } from "react";
import Grid from "@material-ui/core/Grid";

import { CreateName, Rule } from "../../models";
import NameComposer from "./name-composer";


// move to separate module bc global toolbar uses it
const createReadyStatuses = (names: CreateName[]): ReactNode => {
    const counts = { ready: 1, notReady: 1 }
    return (
        <div className="names-status">
            <div>{counts.ready}</div>
            <div>{counts.notReady}</div>
        </div>
    )
}


const NameComposeGroup = ({ names, rule }: { names: CreateName[], rule: Rule }) => {
    return (
        <div className="name-compose-group">
            <Grid container>
                <Grid item xs={3}>
                    <div className="name-compose-group-tools">tools</div>
                </Grid>
                <Grid item xs={3}>
                    <div className="name-compose-group-name">{rule.name}</div>
                </Grid>
                <Grid item xs={3}>
                    {createReadyStatuses(names)}
                </Grid>
            </Grid>
            <ul className="name-compose-list">
                {
                    names.map((name) => {
                        return (
                            <li key={name.localId}>
                                <NameComposer name={name} rule={rule} />
                            </li>
                        )
                    })
                }
            </ul>
        </div>
    )
}

export default NameComposeGroup;