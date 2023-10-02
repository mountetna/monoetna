import React, { useEffect } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector, useDispatch } from 'react-redux'

import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "../names-toolbar/names-toolbar";
import { CreateName } from "../../models";
import CreateNameGroupCompose from "./names-compose-group";
import { fetchRulesFromMagma } from "../../actions/rules";
import { selectCreateNameGroupIdsWithPrimaryRule } from "../../selectors/names";



const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const dispatch = useDispatch();
    const createNameGroupsIdsByPrimaryRule: Record<string, string[]> = useSelector(selectCreateNameGroupIdsWithPrimaryRule)

    const iterateOnSelection = (names: CreateName[], tokenValue: string, start: number, finish: number) => {
        console.log(`adding from ${start} to ${finish} on token "${tokenValue}" with names ${names}`)
    }

    useEffect(() => {
        dispatch(fetchRulesFromMagma(project_name))
    }, []);


    return (
        <>
            <ProjectHeader project_name={project_name} />
            {/* TODO?: render toolbar buttons directly so you don't have to pass props through */}
            <NamesToolbar />
            <Grid className="create-name-group-composers" container>
                {/* TODO: move grouping to selector in redux so we can memoize */}
                {
                    Object.entries(createNameGroupsIdsByPrimaryRule).map(([ruleName, createNameGroupIds]) => {
                        return (
                            <Grid item key={ruleName} xs={4}>
                                <CreateNameGroupCompose
                                    createNameGroupIds={createNameGroupIds}
                                    ruleName={ruleName}
                                />
                            </Grid>
                        )
                    })
                }
            </Grid>
        </>
    )
};

export default NamesCreate;