import React from "react";
import ProjectHeader from "etna-js/components/project-header";
import Grid from "@material-ui/core/Grid"
import Link from "@material-ui/core/Link"
// TODO: use etna navigation
// import Link from "etna-js/components/link"


const ProjectDetail = ({ project_name }: { project_name: string }) => {
    return (
        <>
            <ProjectHeader project_name={project_name} />
            <Grid container>
                <Grid item>
                    <Link href={`/${project_name}/browse`}>Browse</Link>
                </Grid>
                <Grid item>
                    <Link href={`/${project_name}/create`}>Create</Link>
                </Grid>
                <Grid item>
                    <Link href={`/${project_name}/rules`}>Rules</Link>
                </Grid>
            </Grid>
        </>
    )
}

export default ProjectDetail;