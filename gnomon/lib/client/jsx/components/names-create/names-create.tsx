import React, { useEffect } from "react"
import Grid from "@material-ui/core/Grid"
import { makeStyles } from '@material-ui/core/styles';
import { useSelector } from 'react-redux'
import _ from "lodash"

import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "../names-toolbar/names-toolbar";
import CreateNameGroupCompose from "./name-group-composer";
import { fetchRulesFromMagma } from "../../actions/rules";
import { selectCreateNameGroupIdsByPrimaryRule, selectCreateNameGroupsByLocalId, selectFilterCreateNameGroupIds, selectFilterEnabledStatus, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId, selectReplaceCreateNameGroupIds, selectSearchCreateNameGroupIds, selectSelectedCreateNameGroupIds } from "../../selectors/names";
import { useDispatch } from "../../utils/redux";
import Counts from "./counts";
import { clearCreateNameGroupsFilter, clearCreateNameGroupsSelection } from "../../actions/names";



const useStyles = makeStyles((theme) => ({
    projectAndToolbarContainer: {
        display: "flex",
        borderBottom: "1px solid #ccc",
        "& > :first-child, & > :last-child": {
            flex: "1"
        },
        "& > .placeholder": {
            visibility: "hidden",
        },
    },
    countsList: {
        marginTop: "0.5em",
        textAlign: "center"
    },
    readyCounts: {
        display: "inline-block",
        '& .count, & .separator': {
            display: "inline-block",
        },
        "& .separator": {
            margin: "0 0.5em"
        },
        "& .count-not-ready": {
            color: "red",
        },
    },
    selectedAndFilteredCounts: {
        display: "inline-block",
        marginLeft: "4em",
        '& .count, & .separator': {
            display: "inline-block",
        },
        "& .separator": {
            margin: "0 0.5em"
        },
        "& .count-filtered": {
            color: "purple",
        },
        "& .count-selected:hover, .count-filtered:hover": {
            cursor: "pointer",
            textDecoration: "line-through",
        },
    },
    composerList: {
        padding: "0 5em",
        margin: "5.5em 0"
    },
    composer: {
        textAlign: "center",
        padding: "0 2em",
    },
}));


function createCountsList(
    total: number,
    notReady: number,
    selected: number,
    filtered: number,
    classes: Record<string, any>,
    handleClickSelected: () => void,
    handleClickFiltered: () => void,
) {
    const readyCounts = [
        {
            name: "total",
            description: `name${total > 1 ? "s" : ""}`,
            value: total,
        },
        {
            name: "not-ready",
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
            name: "filtered",
            description: "hidden by filter",
            value: filtered,
            hideAtZero: true,
        },
    ]

    return (
        <div className={classes.countsList}>
            <React.Fragment
                key="ready"
            >
                <Counts
                    counts={readyCounts}
                    className={classes.readyCounts}
                    separator="•"
                />
            </React.Fragment>
            {
                _.some(selectedAndFilteredCounts, (count) => count.value > 0)
                    ? (<React.Fragment
                        key="selected-and-filtered"
                    >
                        <Counts
                            counts={selectedAndFilteredCounts}
                            className={classes.selectedAndFilteredCounts}
                            separator="•"
                            onClick={(countName => countName == "selected" ? handleClickSelected() : handleClickFiltered())}
                        />
                    </React.Fragment>)
                    : undefined
            }
        </div>
    )
}


const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const dispatch = useDispatch()
    const classes = useStyles()

    const createNameGroupsIdsByPrimaryRule = useSelector(state => selectCreateNameGroupIdsByPrimaryRule(state, true, true))
    const completeCreateNameGroupsCount = Object.keys(useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)).length
    let totalCreateNameGroupsCount = Object.keys(useSelector(selectCreateNameGroupsByLocalId)).length
    const searchCreateNameGroupIds = useSelector(selectSearchCreateNameGroupIds).size
    const replaceCreateNameGroupIds = useSelector(selectReplaceCreateNameGroupIds).size
    totalCreateNameGroupsCount -= (searchCreateNameGroupIds + replaceCreateNameGroupIds)
    const selectionCreateNameGroupsCount = useSelector(selectSelectedCreateNameGroupIds).size
    const filterCreateNameGroupsCount = useSelector(selectFilterCreateNameGroupIds).size
    const filterEnabled = useSelector(selectFilterEnabledStatus)

    useEffect(() => {
        dispatch(fetchRulesFromMagma(project_name))
    }, []);

    const handleClickSelected = () => {
        dispatch(clearCreateNameGroupsSelection())
    }

    const handleClickFiltered = () => {
        dispatch(clearCreateNameGroupsFilter())
    }

    return (
        <React.Fragment>
            <div className={classes.projectAndToolbarContainer}>
                <ProjectHeader project_name={project_name} />
                <NamesToolbar />
                <ProjectHeader project_name={project_name} className="placeholder" />
            </div>
            {
                totalCreateNameGroupsCount > 0
                    ? createCountsList(
                        totalCreateNameGroupsCount,
                        totalCreateNameGroupsCount - completeCreateNameGroupsCount,
                        selectionCreateNameGroupsCount,
                        filterEnabled ? totalCreateNameGroupsCount - filterCreateNameGroupsCount : 0,
                        classes,
                        handleClickSelected,
                        handleClickFiltered,
                    )
                    : undefined
            }
            <Grid className={classes.composerList} container>
                {
                    Object.entries(createNameGroupsIdsByPrimaryRule).map(([ruleName, createNameGroupLocalIds]) => {
                        return (
                            <Grid item key={ruleName} xs={12} md={6} lg={6} xl={3} className={classes.composer}>
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