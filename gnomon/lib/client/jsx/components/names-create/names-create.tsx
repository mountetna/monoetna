import React, { useState, useEffect } from "react";
import Grid from "@material-ui/core/Grid";
import { v4 as uuidv4 } from 'uuid';
import { useSelector, useDispatch } from 'react-redux'

import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "../names-toolbar/names-toolbar";
import { CreateName } from "../../models";
import NameComposeGroup from "./names-compose-group";
import { selectRules } from "../../selectors/rules";
import { fetchRulesFromMagma } from "../../actions/rules";


interface RuleNameGroupDictionary {
    [index: string]: CreateName[]
}


// TODO: move this to selector in redux
const groupNamesByRule = (names: CreateName[]): RuleNameGroupDictionary => {
    const ruleNameGroups = {} as RuleNameGroupDictionary;

    names.forEach(name => {
        if (!(name.ruleName in ruleNameGroups)) {
            ruleNameGroups[name.ruleName] = []
        }

        ruleNameGroups[name.ruleName].push(name);
    })

    return ruleNameGroups;
}


const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const dispatch = useDispatch();
    const rules = useSelector(selectRules);
    const [names, setNames] = useState<CreateName[]>([]);

    const addNameForRule = (ruleName: string) => {
        console.log(`adding ${ruleName}`)
        const newName = {
            type: "CreateName",
            localId: uuidv4(),
            elements: [],
            ruleName: ruleName,
            selected: false,
        } as CreateName
        setNames((names) => names.concat(newName))
    }

    const iterateOnSelection = (names: CreateName[], tokenValue: string, start: number, finish: number) => {
        console.log(`adding from ${start} to ${finish} on token "${tokenValue}" with names ${names}`)
    }

    useEffect(() => {
        dispatch(fetchRulesFromMagma(project_name))
    }, []);

    const renderComponent = () => {
        return (
            <>
                <ProjectHeader project_name={project_name} />
                {/* TODO?: render toolbar buttons directly so you don't have to pass props through */}
                <NamesToolbar
                    names={names}
                    rules={rules}
                    handleAddNameForRule={addNameForRule}
                    handleAddFromSelection={iterateOnSelection}
                />
                <Grid className="name-compose-groups" container>
                    {/* TODO: move grouping to selector in redux so we can memoize */}
                    {
                        Object.entries(groupNamesByRule(names)).map(([ruleName, nameGroup]) => {
                            const rule = rules.find((rule) => rule.name == ruleName);
                            return (
                                <Grid item key={rule.name} xs={4}>
                                    <NameComposeGroup names={nameGroup} rule={rule} />
                                </Grid>
                            )
                        })
                    }
                </Grid>
            </>
        )
    }

    return (
        <>
            {rules.length ? renderComponent() : "Loading..."}
        </>
    )
};

export default NamesCreate;