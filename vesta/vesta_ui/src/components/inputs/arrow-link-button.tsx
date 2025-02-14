import * as React from 'react';
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import Box from '@mui/material/Box';
import Image from 'next/image';

import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg'


export default function ArrowLinkButton({
    onClick,
    text,
}: {
    onClick: () => void,
    text?: string,
}) {
    const hasText = text !== undefined

    return (
        <ButtonBase
            onClick={onClick}
            sx={{
                display: 'inline-flex',
                alignItems: 'center',
                borderRadius: hasText ? '60px' : '50%',
                bgcolor: 'ground.grade10',
            }}
        >
            {hasText && (
                <Typography
                    variant='pLarge'
                    sx={{
                        display: 'flex',
                        p: '10px',
                        pl: hasText ? '24px' : '10px',
                        color: 'utilityHighlight.main'
                    }}
                >
                    {text}
                </Typography>
            )}

            <Box
                sx={{
                    display: 'flex',
                    p: '10px',
                }}
            >
                <Image
                    src={arrowUpRightLight}
                    alt='Arrow pointing up-right'
                />
            </Box>
        </ButtonBase>
    )
}