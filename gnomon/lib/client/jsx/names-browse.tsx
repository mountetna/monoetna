import React from "react";
import ProjectHeader from "etna-js/components/project-header";

const NamesBrowse = ({ project_name }: { project_name: string }) => {
    return (
        <>
            <ProjectHeader project_name={project_name} />
        </>
    )
};

export default NamesBrowse;