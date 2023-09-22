import React, { useState, useEffect, FunctionComponent } from "react";

import { json_get } from "etna-js/utils/fetch";
import { magmaPath } from "etna-js/api/magma_api";
import ProjectHeader from "etna-js/components/project-header";

import NamesToolbar from "../names-toolbar/names-toolbar";
import { CreateName, Rule, Token, TokenValue, Synonym, Counter, RuleRef } from "../models";
import NameComposeGroup from "./names-compose-group";


interface RuleNameGroupDictionary {
    [index: string]: CreateName[]
}


// TODO: move this to selector in redux
const groupNamesByRule = (names: CreateName[]): RuleNameGroupDictionary => {
    const ruleNameGroups = {} as RuleNameGroupDictionary;

    names.forEach(name => {
        if (!(name.rule_name in ruleNameGroups)) {
            ruleNameGroups[name.rule_name] = []
        }

        ruleNameGroups[name.rule_name].push(name);
    })

    return ruleNameGroups;
}


interface MagmaRulesResponse {
    rules: object
    synonyms: Synonym[]
    tokens: object
}

// move to util module, or make Rule a class
const parseMagmaRulesResponse = (res: MagmaRulesResponse): Rule[] => {
    const tokens = Object.values(res.tokens).map((res_token) => {
        const values = Object.entries(res_token.values).map(([name, label]) => {
            return { name, label } as TokenValue;
        })
        return { name: res_token.name, label: res_token.label, values } as Token;
    });

    const rules: Rule[] = Object.entries(res.rules).map(([name, element_str]) => {
        const elements = element_str.split(" ").map((el: string) => {
            if (el[0] == ".") {
                if (el[1] == "n") { return {} as Counter }

                return { value: el.slice(1) } as RuleRef;
            }
            return tokens.find(token => token.name == el);
        });
        return { name, elements } as Rule;
    });

    return rules;
}


const NamesCreate = ({ project_name }: { project_name: string }) => {
    // TODO: loading state
    const [rules, setRules] = useState<Rule[]>([]);
    const [names, setNames] = useState<CreateName[]>([]);

    const addNameForRule = (rule_name: string) => {
        console.log(`adding ${rule_name}`)
        const newName = { elements: [], rule_name, selected: false } as CreateName
        setNames((names) => names.concat(newName))
    }

    const iterateOnSelection = (names: CreateName[], token_value: string, start: number, finish: number) => {
        console.log(`adding from ${start} to ${finish} on token "${token_value}" with names ${names}`)
    }

    useEffect(() => {
        // TODO: get async/await working
        json_get(magmaPath(`gnomon/${project_name}`)).then(
            ({ config }: { config: MagmaRulesResponse }) => {
                setRules(parseMagmaRulesResponse(config))
            }
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
            <div className="name-compose-groups">
                {/* TODO: move this to selector in redux */}
                {
                    Object.entries(groupNamesByRule(names)).map(([rule_name, nameGroup]) => {
                        const rule = rules.find((rule) => rule.name == rule_name);
                        return <NameComposeGroup names={nameGroup} rule={rule} key={rule.name} />
                    })
                }
            </div>
        </>
    )
};

export default NamesCreate;