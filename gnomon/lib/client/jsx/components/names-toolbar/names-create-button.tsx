import React, { useState } from 'react';
import { useSelector } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Switch from '@material-ui/core/Switch';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import SaveOutlinedIcon from '@material-ui/icons/SaveOutlined';
import AutorenewIcon from '@material-ui/icons/Autorenew';
import CheckCircleOutlineIcon from '@material-ui/icons/CheckCircleOutline';
import ErrorOutlineIcon from '@material-ui/icons/ErrorOutline';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import { useTheme , Theme} from '@material-ui/core/styles';
import _ from 'lodash';

import { useDispatch } from '../../utils/redux';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';
import { selectCompleteCreateNamesCreationPayloads, selectComposeErrorCount, selectNamesCreationRequestState, selectRenderedCompleteCreateNamesByCreateNameGroupLocalId } from '../../selectors/names';
import { makeCreateNamesCreationRequest } from '../../actions/names';
import { selectPathParts } from '../../selectors/location';
import { exportDataToBlob, FILE_FORMATS_TO_MIME } from '../../utils/export';
import { Status } from '../../utils/models';
import MultiOptionButton from '../multi-option-button';
import { isSuperset } from '../../utils/set';
import Table from '../table';



const useStyles = makeStyles((theme) => ({
    dialogRoot: {
        '& [aria-labelledby="create-all-dialog-title"]': {
            padding: '0.5em',
        },
    },
    dialogTitle: {
        paddingBottom: '0',
    },
    dialogContent: {
        '& > *': {
            marginBottom: '1em',
        }
    },
    tableControls: {
        textAlign: 'right',
    },
    showImplicitSwitchLabel: {
        margin: '0',
    },
    table: {
        minHeight: '400px',
        height: '25vh',
        maxHeight: '100%',
        minWidth: '400px',
        width: '25vw',
        maxWidth: '100%',
        [theme.breakpoints.down('sm')]: {
            width: '100%',
            height: '90%',
        },
        '& .explicit': {
            fontWeight: 'bold',
        },
        '& .implicit': {
            fontStyle: 'italic',
            fontWeight: 'normal',
            opacity: '0.5',
        },
    },
    dialogActions: {
        '&.withStatus': {
            justifyContent: 'space-between',
        },
    },
    requestStatus: {
        display: 'inline-flex',
        alignItems: 'center',
        marginRight: '2em',
        '& svg': {
            marginRight: '0.25em',
        },
        '&.inProgress': { color: 'orange', },
        '&.success': { color: 'green', },
        '&.error': { color: 'red', },
    },
    buttonsContainer: {
        '& > button:not(:last-child)': {
            marginRight: '0.75em',
        },
    },
    createButton: {
        '&:not(:disabled)': {
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
    exportButton: {
        '&, & *': {
            fontSize: '13px',
        },
        '& button': {
            padding: '4px 10px',
            minWidth: 'unset',
            '&.menuToggle': {
                padding: '4px 8px',
            },
        },
    },
}));


const NamesCreateButton = ({ small, className }: { small: boolean, className?: string }) => {
    const classes = useStyles();
    const dispatch = useDispatch();
    const theme = useTheme();
    const fullScreen = useMediaQuery(theme.breakpoints.down('sm'));
    const [open, setOpen] = useState<boolean>(false);
    const [showImplicit, setShowImplicit] = useState<boolean>(false);
    const [exportStatus, setExportStatus] = useState<Status>('idle');

    const projectName = useSelector(selectPathParts)[0];
    const completeCreateNameGroupsCount = Object.keys(useSelector(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId)).length;
    const creationRequestPayloads = useSelector(selectCompleteCreateNamesCreationPayloads);
    const creationRequestState = useSelector(selectNamesCreationRequestState);
    const composeErrorCount = useSelector(selectComposeErrorCount);

    let foundImplicit = false;
    const rows = creationRequestPayloads
        .filter(payload => {
            if (payload.implicit) {
                foundImplicit = true;
            }
            return !payload.implicit || showImplicit;
        })
        .map(payload => ({
            name: payload.renderedName,
            rule: payload.ruleName,
            implicit: payload.implicit
        }));

    const handleToggle = () => {
        setOpen(prev => !prev);
    };

    const handleClose = () => {
        // lock out the main view if request made
        if (creationRequestState.status != 'idle') {
            return;
        }
        setOpen(false);
    };

    const handleChangeShowImplicit = (event: React.ChangeEvent<HTMLInputElement>) => {
        setShowImplicit(event.target.checked);
    };

    const handleCreateAll = () => {
        dispatch(
            makeCreateNamesCreationRequest(
                projectName,
                creationRequestPayloads.filter(payload => !payload.implicit)
            )
        );
    };

    const handleExportFile = async (fileFormat: keyof typeof FILE_FORMATS_TO_MIME) => {
        setExportStatus('inProgress');

        try {
            const blob = await exportDataToBlob(
                rows.map(row => ({
                    name: row.name,
                    rule: row.rule,
                })),
                fileFormat,
            );
            const filename = `names-${Date.now()}.${fileFormat}`;

            const a = document.createElement('a');
            a.setAttribute('href', window.URL.createObjectURL(blob));
            a.setAttribute('download', filename);
            a.click();
            a.remove();

            setExportStatus('idle');
        } catch (err) {
            console.error(`Error export names file: ${err}`);
            setExportStatus('error');
        }
    };

    const handleNavigateBack = () => {
        window.location.href = `/${projectName}`;
    };

    const renderCreationRequestStatus = (status: Status) => {
        let icon: JSX.Element | undefined = undefined;

        switch (status) {
            case 'idle':
                return;
            case 'inProgress':
                icon = <AutorenewIcon />;
                break;
            case 'success':
                icon = <CheckCircleOutlineIcon />;
                break;
            case 'error':
                icon = <ErrorOutlineIcon />;
                break;
        }

        return <span className={`${classes.requestStatus} ${creationRequestState.status}`}>
            {icon} {_.startCase(creationRequestState.status)}
        </span>;
    };

    const allNamesCreated = () => {
        if (creationRequestState.status != 'success') {
            return false;
        }

        const expected = new Set(creationRequestPayloads.filter(val => !val.implicit).map(val => val.renderedName));
        const created = new Set((creationRequestState.response?.created || []).map(val => val.name));

        return isSuperset(created, expected);
    };

    return (
        <React.Fragment>
            <ToolbarButtonWithPopper
                text="Create All"
                iconComponent={<SaveOutlinedIcon />}
                variant={small ? 'compact' : 'full'}
                color="#0057FF"
                popperId="create-all-dialog"
                disabled={completeCreateNameGroupsCount - composeErrorCount == 0}
                onClickOrPopperChange={handleToggle}
                className={className}
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
                    {'Create All Names'}
                </DialogTitle>
                <DialogContent className={classes.dialogContent}>
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
                    <Table
                        rows={rows}
                        columns={['name', 'rule']}
                        className={classes.table}
                        getRowClass={row => row.implicit ? 'implicit' : 'explicit'}
                    />
                </DialogContent>
                <DialogActions
                    className={`${classes.dialogActions} ${creationRequestState.status != 'idle' ? 'withStatus' : ''}`}
                >
                    {renderCreationRequestStatus(exportStatus != 'idle' ? exportStatus : creationRequestState.status)}
                    <div className={classes.buttonsContainer}>
                        {
                            creationRequestState.status != 'success' || !allNamesCreated()
                                ? (<React.Fragment>
                                    <Button
                                        onClick={handleClose}
                                        color="secondary"
                                        disableElevation
                                        disabled={creationRequestState.status != 'idle'}
                                        autoFocus
                                    >
                                        Cancel
                                    </Button>
                                    <Button
                                        onClick={handleCreateAll}
                                        color="primary"
                                        disableElevation
                                        disabled={creationRequestState.status != 'idle'}
                                        className={classes.createButton}
                                        autoFocus
                                    >
                                        Create
                                    </Button>
                                </React.Fragment>)
                                : (<React.Fragment>
                                    <Button
                                        onClick={handleNavigateBack}
                                        color="secondary"
                                        disableElevation
                                        autoFocus
                                    >
                                        Back
                                    </Button>
                                    <MultiOptionButton
                                        onClick={handleExportFile}
                                        options={['xlsx', 'csv', 'tsv']}
                                        optionPrefix="Export "
                                        color="primary"
                                        disableElevation
                                        disableRipple
                                        disabled={exportStatus != 'idle'}
                                        className={classes.exportButton}
                                        autoFocus
                                    >
                                        Export
                                    </MultiOptionButton>
                                </React.Fragment>)
                        }
                    </div>
                </DialogActions>
            </Dialog>
        </React.Fragment >
    );
};

export default NamesCreateButton;