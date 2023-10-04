import React, { useState } from "react";
import Grid from "@material-ui/core/Grid";
import { useSelector, useDispatch } from 'react-redux'
import Button from "@material-ui/core/Button";
import FindReplaceOutlinedIcon from "@material-ui/icons/FindReplaceOutlined";
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import SaveOutlinedIcon from "@material-ui/icons/SaveOutlined";
import FilterListIcon from '@material-ui/icons/FilterList';

import AddNamesButton from "./add-names-button";
import AddFromSelectionButton from "./add-from-selection-button";
import { deleteSelectedGroupsWithNames } from "../../actions/names";
import { selectSelectedCreateNameGroupIds } from "../../selectors/names";



const NamesToolbar = () => {
    const dispatch = useDispatch()
    const [small, setSmall] = useState<boolean>(false);

    const selectedCreateNameGroupLocalIds: string[] = useSelector(selectSelectedCreateNameGroupIds)

    const handleClickDeleteSelection = () => {
        dispatch(deleteSelectedGroupsWithNames())
    }

    return (
        <Grid container>
            <Grid item xs={2}>
                <AddNamesButton />
            </Grid>
            <Grid item xs={2}>
                <AddFromSelectionButton />
            </Grid>
            <Grid item xs={2}>
                <Button
                    startIcon={<FindReplaceOutlinedIcon />}
                    color="primary"
                    disableElevation
                    disabled={selectedCreateNameGroupLocalIds.length == 0}
                >
                    Replace in Selection
                </Button>
            </Grid>
            <Grid item xs={2}>
                <Button
                    onClick={handleClickDeleteSelection}
                    startIcon={<DeleteOutlineOutlinedIcon />}
                    color="primary"
                    disableElevation
                    disabled={selectedCreateNameGroupLocalIds.length == 0}
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