import React, { useState } from "react";
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";

import AddNamesButton from "./add-names-button";
import FindAndFilterButton from "./find-and-filter-button";
import CopyAndReplaceButton from "./copy-and-replace-button"
import NamesCreateButton from "./names-create-button";
import DeleteButton from "./delete-button";
import GrowShrinkButton from "./grow-shrink-button";



const useStyles = makeStyles((theme) => ({
    container: {
        textAlign: "center",
    },
    buttonsOuterContainer: {
        display: "inline-block",
        // border: "1px solid #ccc",
        borderTop: "none",
    },
    buttonsContainer: {
        display: "inline-block",
        padding: "1.25em 0",
    },
    growShrinkButton: {
        "& svg": {
            fontSize: "1rem",
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

                    <GrowShrinkButton
                        small={small}
                        onClick={() => setSmall(prev => !prev)}
                    />
                </div>
            </div>
        </div>
    )
};

export default NamesToolbar;