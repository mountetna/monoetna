import React, { useState } from "react";
import { makeStyles } from '@material-ui/core/styles';

import GrowShrinkButton from "../names-toolbar/grow-shrink-button";



const useStyles = makeStyles((theme) => ({
    container: {
        textAlign: "center",
    },
    buttonsOuterContainer: {
        display: "inline-block",
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


interface ButtonProps {
    small: boolean
}


const NamesToolbar = ({ buttons, canBeSmall = true }: {
    buttons: React.ReactElement<ButtonProps>[],
    canBeSmall?: boolean,
}) => {
    const classes = useStyles()

    const [small, setSmall] = useState<boolean>(canBeSmall);

    return (
        <div className={classes.container}>
            <div className={classes.buttonsOuterContainer}>
                <div className={classes.buttonsContainer}>
                    {buttons.map((ButtonEl, idx) =>
                        <React.Fragment key={idx}>
                            {React.cloneElement(ButtonEl, { small })}
                        </React.Fragment>
                    )}
                    <GrowShrinkButton
                        small={small}
                        onClick={() => canBeSmall && setSmall(prev => !prev)}
                    />
                </div>
            </div>
        </div>
    )
};

export default NamesToolbar;