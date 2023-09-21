import React, { useState, useEffect } from "react";

import { json_get } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';
import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "./names-toolbar/toolbar";
import { CreateName, SelectableToken } from "./models";



const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const [rules, setRules] = useState<string[]>([]);
    const [names, setNames] = useState<CreateName[]>([]);

    const addNameForRule = (rule: string) => {
        console.log(`adding ${rule}`)
    }

    const iterateOnSelection = (names: CreateName[], token: SelectableToken, start: number, finish: number) => {
        console.log(`adding from ${start} to ${finish} on token "${token}" with names ${names}`)
    }

    useEffect(() => {
        // TODO: get async/await working
        json_get(magmaPath(`gnomon/${project_name}`)).then(
            ({ config }) => setRules(config.rules)
        )
    }, []);

    return (
        <>
            <ProjectHeader project_name={project_name} />
            <NamesToolbar
                names={names}
                rules={rules}
                handleAddNameForRule={addNameForRule}
                handleAddFromSelection={iterateOnSelection}
            />
        </>
    )
};

export default NamesCreate;