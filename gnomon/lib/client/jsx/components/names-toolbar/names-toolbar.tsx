import React, { useState, useRef } from "react";
import Grid from "@material-ui/core/Grid";
import Button from "@material-ui/core/Button";
import FindReplaceOutlinedIcon from "@material-ui/icons/FindReplaceOutlined";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import SaveOutlinedIcon from "@material-ui/icons/SaveOutlined";
import FilterListIcon from '@material-ui/icons/FilterList';

import { CreateName, Rule } from "../../models";
import AddNamesButton from "./add-names-button";
import AddFromSelectionButton from "./add-from-selection-button";



const NamesToolbar = ({ rules, names, handleAddNameForRule, handleAddFromSelection }:
    {
        names: CreateName[],
        rules: Rule[],
        handleAddNameForRule: (rule_name: string) => any,
        handleAddFromSelection: (names: CreateName[], tokenValue: string, start: number, finish: number) => any
    }) => {

    const [small, setSmall] = useState<Boolean>(false);

    return (
        <Grid container>
            <Grid item xs={2}>
                <AddNamesButton
                    ruleNames={rules.map((rule) => rule.name)}
                    clickAddHandler={handleAddNameForRule}
                />
            </Grid>
            <Grid item xs={2}>
                <AddFromSelectionButton
                    selection={names.filter(name => name.selected)}
                    rules={rules}
                    clickAddHandler={handleAddFromSelection}
                />
            </Grid>
            <Grid item xs={2}>
                <Button
                    startIcon={<FindReplaceOutlinedIcon />}
                    color="primary"
                    disableElevation
                >
                    Replace in Selection
                </Button>
            </Grid>
            <Grid item xs={2}>
                <Button
                    startIcon={<DeleteOutlineOutlinedIcon />}
                    color="primary"
                    disableElevation
                >
                    Delete Selection
                </Button>
            </Grid>
            <Grid item xs={2}>
                <Button
                    startIcon={<FilterListIcon />}
                    color="primary"
                    disableElevation
                >
                    Filter
                </Button>
            </Grid>
            <Grid item xs={2}>
                <Button
                    startIcon={<SaveOutlinedIcon />}
                    color="secondary"
                    disableElevation
                >
                    Create All
                </Button>
            </Grid>
        </Grid>
    )
};

export default NamesToolbar;