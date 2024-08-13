'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { SxProps, useTheme } from '@mui/material';
import Image from 'next/image';

import { PrincipalInvestigator } from './models';


// Maybe make this a generic member component
// instead of specific to PIs?

// TODO: add profile link if available
export default function ProjectPI({
    data,
    showAvatar = true,
    showNameAndTitle = true,
    variant = 'stroked',
}: {
    data: PrincipalInvestigator,
    showAvatar?: boolean,
    showNameAndTitle?: boolean,
    variant?: 'stroked' | 'filled',
}) {
    const theme = useTheme()

    const [imageLoaded, setImageLoaded] = React.useState(false)

    let rootStyles: SxProps = {}
    if (showNameAndTitle) {
        rootStyles = {
            p: '4px',
            pr: '18px',
            ...(variant === 'filled' ? {} : {
                bgcolor: 'utilityWhite.main',
                border: `1px solid ${theme.palette.ground.grade75}`,
            })
        }
    }

    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'row',
                gap: '5.13px',
                borderRadius: '70px',
                ...rootStyles,
            }}
        >
            <Box
                sx={{
                    minWidth: '41px',
                    maxWidth: '41px',
                    minHeight: '41px',
                    maxHeight: '41px',
                    borderRadius: '50%',
                    overflow: 'hidden',
                    position: 'relative',
                    '& > *': {
                        width: '100%',
                        minWidth: '100%',
                        height: '100%',
                        minHeight: '100%',
                    },
                }}
            >
                <Box
                    sx={{
                        display: 'flex',
                        justifyContent: 'center',
                        alignItems: 'center',
                        bgcolor: data.color,
                    }}
                >
                    <Typography
                        variant='pBodyBoldWt'
                        sx={{
                            color: data.altColor,
                        }}
                    >
                        {data.name[0].toUpperCase()}
                    </Typography>
                </Box>

                {showAvatar && data.imageUrl &&
                    <Image
                        src={data.imageUrl}
                        alt={`Profile image for ${data.name}`}
                        width={41}
                        height={41}
                        onLoad={() => setImageLoaded(true)}
                        style={{
                            opacity: imageLoaded ? 1 : 0,
                            position: 'absolute',
                            top: 0,
                            left: 0,
                            transition: theme.transitions.create(
                                ['opacity'],
                                {
                                    easing: theme.transitions.easing.swell,
                                    duration: theme.transitions.duration.swell,
                                },
                            ),
                        }}
                    />
                }
            </Box>

            {showNameAndTitle && <Box
                sx={{
                    overflow: 'hidden',
                    // justifyContent: 'space-around',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center',
                    '& *': {
                        display: 'block',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        // textWrap: 'nowrap',
                        textWrap: 'wrap',
                        // '&:hover, &:focus': {
                        //     textWrap: 'nowrap',
                        // },
                    },
                }}
            >
                <Typography variant='pBodyBoldWt'>
                    {data.name}
                </Typography>

                {data.title &&
                    <Typography variant='p3XS'>
                        {data.title}
                    </Typography>
                }
            </Box>}
        </Box>
    )
}