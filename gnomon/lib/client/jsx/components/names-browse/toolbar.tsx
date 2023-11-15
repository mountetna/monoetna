import React, { useState } from "react";
import { makeStyles } from '@material-ui/core/styles';

import GrowShrinkButton from "../names-toolbar/grow-shrink-button";
import ExportButton from "../names-toolbar/export-button";



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


const NamesBrowseToolbar = ({ exportData, exportButtonText }: { exportData: Array<any>, exportButtonText: string }) => {
    const classes = useStyles()

    const [small, setSmall] = useState<boolean>(true);

    return (
        <div className={classes.container}>
            <div className={classes.buttonsOuterContainer}>
                <div className={classes.buttonsContainer}>
                    <ExportButton
                        small={small}
                        data={exportData}
                        buttonText={exportButtonText}
                    />
                    <GrowShrinkButton
                        small={small}
                        onClick={() => setSmall(prev => !prev)}
                    />
                </div>
            </div>
        </div>
    )
};

export default NamesBrowseToolbar;