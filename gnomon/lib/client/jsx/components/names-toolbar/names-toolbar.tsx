import React, { useState } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector } from 'react-redux'
import Button from "@material-ui/core/Button";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import SaveOutlinedIcon from "@material-ui/icons/SaveOutlined";

import AddNamesButton from "./add-names-button";
import { deleteSelectedGroupsWithNames } from "../../actions/names";
import { selectSelectedCreateNameGroupIds } from "../../selectors/names";
import FindAndFilterButton from "./find-and-filter-button";
import ReplaceInSelectionButton from "./copy-and-replace-button"
import { useDispatch } from "../../utils/redux";
import { State } from "../../store";
import { selectGlobalState } from "../../selectors/global";



const NamesToolbar = () => {
    const dispatch = useDispatch()
    const [small, setSmall] = useState<boolean>(false);

    const globalState: State = useSelector(selectGlobalState)
    const selectedCreateNameGroupLocalIds: Set<string> = useSelector(selectSelectedCreateNameGroupIds)

    const handleClickDelete = () => {
        dispatch(deleteSelectedGroupsWithNames(globalState))
    }

    return (
        <Grid container>
            <Grid item xs={2}>
                <AddNamesButton />
            </Grid>
            <Grid item xs={2}>
                <FindAndFilterButton />
            </Grid>
            <Grid item xs={2}>
                <ReplaceInSelectionButton />
            </Grid>
            <Grid item xs={2}>
                <Button
                    onClick={handleClickDelete}
                    startIcon={<DeleteOutlineOutlinedIcon />}
                    color="primary"
                    disableElevation
                    disabled={selectedCreateNameGroupLocalIds.size == 0}
                >
                    Delete
                </Button>
            </Grid>
            <Grid item xs={2}>
                <Button
                    startIcon={<SaveOutlinedIcon />}
                    color="secondary"
                    disableElevation
                    disabled={true}
                >
                    Create All
                </Button>
            </Grid>
        </Grid>
    )
};

export default NamesToolbar;