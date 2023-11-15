import React, { useState } from "react";
import { makeStyles } from '@material-ui/core/styles';

import AddNamesButton from "../names-toolbar/add-names-button";
import FindAndFilterButton from "../names-toolbar/find-and-filter-button";
import CopyAndReplaceButton from "../names-toolbar/copy-and-replace-button"
import NamesCreateButton from "../names-toolbar/names-create-button";
import DeleteButton from "../names-toolbar/delete-button";
import GrowShrinkButton from "../names-toolbar/grow-shrink-button";



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


const NamesCreateToolbar = () => {
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

export default NamesCreateToolbar;