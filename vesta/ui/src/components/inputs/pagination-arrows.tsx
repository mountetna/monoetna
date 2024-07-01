import * as React from 'react'
import Image from "next/image"
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';

import arrowUpRightDarkSrc from '/public/images/icons/arrow-up-right-dark.svg'
import { SxProps, useTheme } from '@mui/material';


export default function PaginationArrows({
    onClickPrev,
    onClickNext,
    prevDisabled = false,
    nextDisabled = false,
}: {
    onClickPrev: () => void,
    onClickNext: () => void,
    prevDisabled?: boolean,
    nextDisabled?: boolean,
}) {
    const theme = useTheme()
    const buttonStyles = (disabled: boolean): SxProps => ({
        color: '',
        bgcolor: disabled ? theme.palette.ground.grade75 : theme.palette.utilityWhite.main,
        transition: theme.transitions.create(
            ['background-color'],
            {
                duration: theme.transitions.duration.ease,
                easing: theme.transitions.easing.ease,
            }
        ),
        p: '8px',
        borderRadius: '50%',
        '& img': {
            width: '24px',
            height: '24px',
        },
    })

    return (
        <Box
        sx={{
            display: 'flex',
            gap: '8px',
        }}
        >
            <ButtonBase
                onClick={onClickPrev}
                disabled={prevDisabled}
                sx={buttonStyles(prevDisabled)}
            >
                <Image
                    src={arrowUpRightDarkSrc}
                    alt='Arrow left'
                    style={{
                        rotate: '225deg',
                    }}
                />
            </ButtonBase>
            <ButtonBase
                onClick={onClickNext}
                disabled={nextDisabled}
                sx={buttonStyles(nextDisabled)}
            >
                <Image
                    src={arrowUpRightDarkSrc}
                    alt='Arrow right'
                    style={{
                        rotate: '45deg',
                    }}
                />
            </ButtonBase>
        </Box>
    )
}