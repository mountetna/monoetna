import React, { useState } from 'react';
import { makeStyles } from '@material-ui/core/styles';

import GrowShrinkButton from '../names-toolbar/grow-shrink-button';



const useStyles = makeStyles((theme) => ({
    container: {
        textAlign: 'center',
    },
    buttonsOuterContainer: {
        display: 'inline-block',
        borderTop: 'none',
    },
    buttonsContainer: {
        display: 'inline-block',
        padding: '1.25em 0',
        [theme.breakpoints.down('sm')]: {
            display: 'flex',
            flexWrap: 'wrap',
            justifyContent: 'center',
        },
    },
    buttonContainer: {
        maxWidth: '8em',
        marginRight: '1.25em',
        '&:last-child': {
            marginRight: '0',
        },
        [theme.breakpoints.down('sm')]: {
            maxWidth: 'unset',
            margin: '0',
            width: 'calc(100% / 3)',
        },
        '&:nth-child(-n+3)': {
            [theme.breakpoints.down('sm')]: {
                paddingBottom: '1.25em',
            },
        },
    },
    growShrinkButton: {
        '& svg': {
            fontSize: '1rem',
        }
    },
}));


interface ButtonProps {
    small: boolean
    className: string
}


const NamesToolbar = ({ buttons, canBeSmall = true }: {
    buttons: React.ReactElement<ButtonProps>[],
    canBeSmall?: boolean,
}) => {
    const classes = useStyles();

    const [small, setSmall] = useState<boolean>(canBeSmall);

    return (
        <div className={classes.container}>
            <div className={classes.buttonsOuterContainer}>
                <div className={classes.buttonsContainer}>
                    {buttons.map((ButtonEl, idx) =>
                        <React.Fragment key={idx}>
                            {React.cloneElement(ButtonEl, { small, className: classes.buttonContainer })}
                        </React.Fragment>
                    )}
                    <GrowShrinkButton
                        small={small}
                        onClick={() => canBeSmall && setSmall(prev => !prev)}
                        className={classes.buttonContainer}
                    />
                </div>
            </div>
        </div>
    );
};

export default NamesToolbar;