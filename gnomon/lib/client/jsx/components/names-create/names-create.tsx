import React, { useEffect } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector } from 'react-redux'

import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "../names-toolbar/names-toolbar";
import CreateNameGroupCompose from "./names-compose-group";
import { fetchRulesFromMagma } from "../../actions/rules";
import { selectCreateNameGroupIdsByPrimaryRule } from "../../selectors/names";
import { useDispatch } from "../../utils/redux";



const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const dispatch = useDispatch();
    const createNameGroupsIdsByPrimaryRule: Record<string, string[]> = useSelector(state => selectCreateNameGroupIdsByPrimaryRule(state, true, true))

    useEffect(() => {
        dispatch(fetchRulesFromMagma(project_name))
    }, []);


    return (
        <>
            <ProjectHeader project_name={project_name} />
            <NamesToolbar />
            <Grid className="create-name-group-composers" container>
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