import React, { useState } from "react";
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import UnfoldLessOutlinedIcon from '@material-ui/icons/UnfoldLessOutlined';
import UnfoldMoreOutlinedIcon from '@material-ui/icons/UnfoldMoreOutlined';
import Tooltip from '@material-ui/core/Tooltip';

import AddNamesButton from "./add-names-button";
import FindAndFilterButton from "./find-and-filter-button";
import CopyAndReplaceButton from "./copy-and-replace-button"
import NamesCreateButton from "./names-create-button";
import DeleteButton from "./delete-button";



const useStyles = makeStyles((theme) => ({
    container: {
        textAlign: "center",
    },
    buttonsOuterContainer: {
        display: "inline-block",
        border: "1px solid #ccc",
        borderTop: "none",
    },
    buttonsContainer: {
        display: "inline-block",
        padding: "1.25em 0",
    },
    growShrinkButton: {
        "& svg": {
            fontSize: "1rem",
            transform: "rotate(90deg)",
        }
    },
}));


const NamesToolbar = () => {
    const classes = useStyles()

    const [small, setSmall] = useState<boolean>(true);

    return (
        <div className={classes.container}>
            <div className={classes.buttonsOuterContainer}>
                <div className={classes.buttonsContainer}>
                    <AddNamesButton small={small} />

                    <FindAndFilterButton small={small} />

                    <CopyAndReplaceButton small={small} />

                    <DeleteButton small={small} />

                    <NamesCreateButton small={small} />
                </div>
                <Button
                    onClick={() => setSmall(prev => !prev)}
                    color="secondary"
                    aria-label={"Toggle Grow/Shrink Toolbar"}
                    disableRipple
                    disableTouchRipple
                    fullWidth={true}
                    size="small"
                    className={classes.growShrinkButton}
                >
                    {small ? <UnfoldMoreOutlinedIcon /> : <UnfoldLessOutlinedIcon />}
                </Button>
            </div>
        </div>
    )
};

export default NamesToolbar;