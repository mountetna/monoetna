import React, { useState } from "react";
import { useSelector } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Switch from '@material-ui/core/Switch';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import SaveOutlinedIcon from "@material-ui/icons/SaveOutlined";
import AutorenewIcon from '@material-ui/icons/Autorenew';
import CheckCircleOutlineIcon from '@material-ui/icons/CheckCircleOutline';
import ErrorOutlineIcon from '@material-ui/icons/ErrorOutline';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import { useTheme } from '@material-ui/core/styles';
import _ from "lodash"
import Excel from "exceljs"

import { useDispatch } from "../../utils/redux";
import NamesTable from "./names-create-table";
import ToolbarButtonWithPopper from "./toolbar-button-with-popper";
import { selectCompleteCreateNamesCreationPayloads, selectNamesCreationRequestState, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId } from "../../selectors/names";
import { makeCreateNamesCreationRequest } from "../../actions/names";
import { selectPathParts } from "../../selectors/location"



const useStyles = makeStyles((theme) => ({
    dialogRoot: {
        '& [aria-labelledby="create-all-dialog-title"]': {
            padding: "0.5em",
        },
    },
    dialogTitle: {
        paddingBottom: "0",
    },
    tableControls: {
        marginBottom: "1em",
        textAlign: "right",
    },
    showImplicitSwitchLabel: {
        margin: "0",
    },
    dialogActions: {
        "&.withStatus": {
            justifyContent: "space-between",
        },
    },
    requestStatus: {
        display: "inline-flex",
        alignItems: "center",
        marginRight: "2em",
        "& svg": {
            marginRight: "0.25em",
        },
        "&.inProgress": { color: "orange", },
        "&.success": { color: "green", },
        "&.error": { color: "red", },
    },
    buttonsContainer: {
        "& > button:not(:last-child)": {
            marginRight: "0.75em",
        },
    },
    createButton: {
        "&:not(:disabled)": {
            backgroundImage: 'url(/images/distort.svg)',
            '&:hover': {
                backgroundImage: 'url(/images/distort2.svg)',
                boxShadow: '0 0 20px #d95'
            },
            '&:active': {
                backgroundImage: 'url(/images/distort3.svg)',
                boxShadow: '0 0 2px #d95'
            },
            backgroundSize: 'cover',
            border: '2px solid darkgoldenrod',
            boxShadow: '0 0 10px #d95',
            color: 'white',
            fontWeight: 'bold'
        },
    },
}));


const NamesCreateButton = ({ small }: { small: boolean }) => {
    const classes = useStyles()
    const dispatch = useDispatch()
    const theme = useTheme();
    const fullScreen = useMediaQuery(theme.breakpoints.down('sm'));
    const [open, setOpen] = useState<boolean>(false);
    const [showImplicit, setShowImplicit] = useState<boolean>(false);

    const projectName = useSelector(selectPathParts)[0]
    const completeCreateNameGroupsCount = Object.keys(useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)).length
    const creationRequestPayloads = useSelector(selectCompleteCreateNamesCreationPayloads)
    const creationRequestState = useSelector(selectNamesCreationRequestState)

    let foundImplicit = false
    const rows = creationRequestPayloads
        .filter(payload => {
            if (payload.implicit) {
                foundImplicit = true
            }
            return !payload.implicit || showImplicit
        })
        .map(payload => ({
            name: payload.renderedName,
            rule: payload.ruleName,
            implicit: payload.implicit
        }))

    const handleToggle = () => {
        setOpen(prev => !prev);
    }

    const handleClose = () => {
        // lock out the main view if request made
        if (creationRequestState.status != "idle") {
            return
        }
        setOpen(false);
    }

    const handleChangeShowImplicit = (event: React.ChangeEvent<HTMLInputElement>) => {
        setShowImplicit(event.target.checked)
    }

    const handleCreateAll = () => {
        dispatch(
            makeCreateNamesCreationRequest(
                projectName,
                creationRequestPayloads.filter(payload => !payload.implicit)
            )
        )
    }

    const renderCreationRequestStatus = () => {
        const status = creationRequestState.status
        let icon: JSX.Element | undefined = undefined

        switch (status) {
            case "idle":
                return
            case "inProgress":
                icon = <AutorenewIcon />
                break
            case "success":
                icon = <CheckCircleOutlineIcon />
                break
            case "error":
                icon = <ErrorOutlineIcon />
                break
        }

        return <span className={`${classes.requestStatus} ${creationRequestState.status}`}>
            {icon} {_.capitalize(creationRequestState.status)}
        </span>
    }

    return (
        <React.Fragment>
            <ToolbarButtonWithPopper
                text="Create All"
                iconComponent={<SaveOutlinedIcon />}
                variant={small ? "compact" : "full"}
                color="#0057FF"
                popperId="create-all-dialog"
                disabled={completeCreateNameGroupsCount == 0}
                onClickOrPopperChange={handleToggle}
            />
            <Dialog
                id="create-all-dialog"
                className={classes.dialogRoot}
                fullScreen={fullScreen}
                open={open}
                onClose={handleClose}
                aria-labelledby="create-all-dialog-title"
            >
                <DialogTitle
                    id="create-all-dialog-title"
                    className={classes.dialogTitle}
                >
                    {"Create All Names"}
                </DialogTitle>
                <DialogContent>
                    <div className={classes.tableControls}>
                        {foundImplicit &&
                            <FormControlLabel
                                className={classes.showImplicitSwitchLabel}
                                control={
                                    <Switch
                                        checked={showImplicit}
                                        onChange={handleChangeShowImplicit}
                                        name="showImplicit"
                                        disableRipple
                                    />
                                }
                                label="Show Implicit"
                            />}
                    </div>
                    <NamesTable
                        rows={rows}
                    />
                </DialogContent>
                <DialogActions
                    className={`${classes.dialogActions} ${creationRequestState.status != "idle" ? "withStatus" : ""}`}
                >
                    {renderCreationRequestStatus()}
                    <div className={classes.buttonsContainer}>
                        {
                            creationRequestState.status == ""
                                ? (<div></div>)
                                : (<React.Fragment>
                                    <Button
                                        autoFocus
                                        onClick={handleClose}
                                        color="secondary"
                                        disableElevation
                                        disabled={creationRequestState.status != "idle"}
                                    >
                                        Cancel
                                    </Button>
                                    <Button
                                        autoFocus
                                        onClick={handleCreateAll}
                                        color="primary"
                                        disableElevation
                                        disabled={creationRequestState.status != "idle"}
                                        className={classes.createButton}
                                    >
                                        Create
                                    </Button>
                                </React.Fragment>)
                        }
                    </div>
                </DialogActions>
            </Dialog>
        </React.Fragment >
    )
}

export default NamesCreateButton