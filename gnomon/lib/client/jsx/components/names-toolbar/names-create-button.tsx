import React, { useState, useRef } from "react";
import { useSelector, batch } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import Switch from '@material-ui/core/Switch';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import SaveOutlinedIcon from "@material-ui/icons/SaveOutlined";
import useMediaQuery from '@material-ui/core/useMediaQuery';
import { useTheme } from '@material-ui/core/styles';

import { useDispatch } from "../../utils/redux";
import NamesTable, { Data } from "../names-table";
import { CompleteCreateNameRequestPayload, selectCompleteCreateNamesCreationPayloads, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId } from "../../selectors/names";



const useStyles = makeStyles((theme) => ({
    container: {
        display: "inline-block",
    },
    tableControls: {
        textAlign: 'right',
    },
}));


const NamesCreateButton = () => {
    const classes = useStyles()
    const dispatch = useDispatch()
    const theme = useTheme();
    const fullScreen = useMediaQuery(theme.breakpoints.down('sm'));
    const [open, setOpen] = useState<boolean>(false);
    const [showImplicit, setShowImplicit] = useState<boolean>(false);

    const completeCreateNameGroupsCount = Object.keys(useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)).length
    const creationRequestPayloads: CompleteCreateNameRequestPayload[] = useSelector(selectCompleteCreateNamesCreationPayloads)
    const rows: Data[] = creationRequestPayloads
        .filter(payload => !payload.implicit || showImplicit)
        .map(payload => { return { name: payload.renderedName, rule: payload.ruleName } })

    const handleToggle = () => {
        setOpen(prev => !prev);
    }

    const handleClose = () => {
        setOpen(false);
    }

    const handleChangeShowImplicit = (event: React.ChangeEvent<HTMLInputElement>) => {
        setShowImplicit(event.target.checked)
    }

    return (
        <div className={classes.container}>
            <Button
                startIcon={<SaveOutlinedIcon />}
                onClick={handleToggle}
                color="secondary"
                aria-label="Open Create All Dialog"
                aria-haspopup="true"
                aria-controls={open ? "create-all-dialog" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
                disabled={completeCreateNameGroupsCount == 0}
            >
                Create All
            </Button>
            <Dialog
                id="create-all-dialog"
                fullScreen={fullScreen}
                open={open}
                onClose={handleClose}
                aria-labelledby="create-all-dialog-title"
            >
                <DialogTitle id="create-all-dialog-title">{"Create All Names"}</DialogTitle>
                <DialogContent>
                    <div className={classes.tableControls}>
                        <FormControlLabel
                            control={
                                <Switch
                                    checked={showImplicit}
                                    onChange={handleChangeShowImplicit}
                                    name="showImplicit"
                                    disableRipple
                                />
                            }
                            label="Show Implicit"
                        />
                    </div>
                    <NamesTable
                        rows={rows}
                    />
                </DialogContent>
                <DialogActions>
                    <Button
                        autoFocus
                        onClick={handleClose}
                        color="secondary"
                        disableElevation
                    >
                        Cancel
                    </Button>
                    <Button
                        autoFocus
                        onClick={handleClose}
                        color="primary"
                        disableElevation
                    >
                        Create All
                    </Button>
                </DialogActions>
            </Dialog>
        </div>
    )
}

export default NamesCreateButton