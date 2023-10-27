import React, { useEffect } from "react"
import Grid from "@material-ui/core/Grid"
import { makeStyles } from '@material-ui/core/styles';
import { useSelector } from 'react-redux'
import * as _ from "lodash"

import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "../names-toolbar/names-toolbar";
import CreateNameGroupCompose from "./names-compose-group";
import { fetchRulesFromMagma } from "../../actions/rules";
import { selectCreateNameGroupIdsByPrimaryRule, selectCreateNameGroupsByLocalId, selectFilterCreateNameGroupIds, selectFilterEnabledStatus, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId, selectReplaceCreateNameGroupIds, selectSearchCreateNameGroupIds, selectSelectedCreateNameGroupIds } from "../../selectors/names";
import { useDispatch } from "../../utils/redux";
import Counts from "./counts";



function createCountsList(
    total: number,
    notReady: number,
    selected: number,
    filtered: number,
    classes: Record<string, any>,
) {
    const readyCounts = [
        {
            name: "total",
            description: `name${total > 1 ? "s" : ""}`,
            value: total,
        },
        {
            name: "notReady",
            description: "not ready",
            value: notReady,
            hideAtZero: true,
        },
    ]
    const selectedAndFilteredCounts = [
        {
            name: "selected",
            description: "selected",
            value: selected,
            hideAtZero: true,
        },
        {
            name: "filter",
            description: "hidden by filter",
            value: filtered,
            hideAtZero: true,
        },
    ]

    return (
        <div className="counts-list">
            <React.Fragment
                key="ready"
            >
                <Counts
                    counts={readyCounts}
                    className={classes.readyCounts}
                />
            </React.Fragment>
            {
                _.some(selectedAndFilteredCounts, (count) => count.value > 0)
                    ? (<React.Fragment
                        key="selection-and-filter"
                    >
                        <Counts
                            counts={selectedAndFilteredCounts}
                            className={classes.selectedAndFilteredCounts}
                        />
                    </React.Fragment>)
                    : undefined
            }
        </div>
    )
}


const useStyles = makeStyles((theme) => ({
    readyCounts: {
        "& count-total": {
            display: "inline-block",
        }
    },
    selectedAndFilteredCounts: {
    },
}));


const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const dispatch = useDispatch()
    const classes = useStyles()

    const createNameGroupsIdsByPrimaryRule: Record<string, string[]> = useSelector(state => selectCreateNameGroupIdsByPrimaryRule(state, true, true))
    const completeCreateNameGroupsCount = Object.keys(useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)).length
    let totalCreateNameGroupsCount = Object.keys(useSelector(selectCreateNameGroupsByLocalId)).length
    const searchCreateNameGroupIds = useSelector(selectSearchCreateNameGroupIds).size
    const replaceCreateNameGroupIds = useSelector(selectReplaceCreateNameGroupIds).size
    totalCreateNameGroupsCount -= (searchCreateNameGroupIds + replaceCreateNameGroupIds)
    const selectionCreateNameGroupsCount = useSelector(selectSelectedCreateNameGroupIds).size
    const filterCreateNameGroupsCount = useSelector(selectFilterCreateNameGroupIds).size
    const filterEnabled: boolean = useSelector(selectFilterEnabledStatus)

    useEffect(() => {
        dispatch(fetchRulesFromMagma(project_name))
    }, []);


    return (
        <React.Fragment>
            <ProjectHeader project_name={project_name} />
            <NamesToolbar />
            {
                totalCreateNameGroupsCount > 0
                    ? createCountsList(
                        totalCreateNameGroupsCount,
                        totalCreateNameGroupsCount - completeCreateNameGroupsCount,
                        selectionCreateNameGroupsCount,
                        filterEnabled ? totalCreateNameGroupsCount - filterCreateNameGroupsCount : 0,
                        classes,
                    )
                    : undefined
            }
            <Grid className="create-name-group-composers" container>
                {
                    Object.entries(createNameGroupsIdsByPrimaryRule).map(([ruleName, createNameGroupLocalIds]) => {
                        return (
                            <Grid item key={ruleName} xs={12} md={6} lg={4} xl={3}>
                                <CreateNameGroupCompose
                                    createNameGroupLocalIds={createNameGroupLocalIds}
                                    ruleName={ruleName}
                                />
                            </Grid>
                        )
                    })
                }
            </Grid>
        </React.Fragment>
    )
};

export default NamesCreate;