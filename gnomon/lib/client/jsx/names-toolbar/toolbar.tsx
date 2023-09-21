import React, { useState, useRef } from "react";
import Grid from "@material-ui/core/Grid";
import Button from "@material-ui/core/Button";

import FindReplaceOutlinedIcon from "@material-ui/icons/FindReplaceOutlined";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import SaveOutlinedIcon from "@material-ui/icons/SaveOutlined";

import { CreateName, SelectableToken, Rule } from "../models";
import AddNamesButton from "./add-names";
import AddFromSelectionButton from "./add-from-selection";



const NamesToolbar = ({ rules, names, handleAddNameForRule, handleAddFromSelection }:
    {
        names: CreateName[],
        rules: Rule[],
        handleAddNameForRule: (rule: string) => any,
        handleAddFromSelection: (names: CreateName[], token: SelectableToken, start: number, finish: number) => any
    }) => {

    const [small, setSmall] = useState<Boolean>(false);

    return (
        <Grid container>
            <Grid item xs={12}>
                <AddNamesButton
                    rules={Object.keys(rules || {})}
                    clickAddHandler={handleAddNameForRule}
                />
                <AddFromSelectionButton
                    selection={names.filter(name => name.selected)}
                    rules={rules}
                    clickAddHandler={handleAddFromSelection}
                />
                <Button
                    startIcon={<FindReplaceOutlinedIcon />}
                    color="primary"
                    disableElevation
                >
                    Replace in Selection
                </Button>
                <Button
                    startIcon={<DeleteOutlineOutlinedIcon />}
                    color="primary"
                    disableElevation
                >
                    Delete Selection
                </Button>
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