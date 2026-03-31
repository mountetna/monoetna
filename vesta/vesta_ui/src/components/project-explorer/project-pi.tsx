'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { SxProps, useTheme } from '@mui/material';
import Image from 'next/image';

import { PrincipalInvestigator } from './models';
import Link from '../link/link';


// Maybe make this a generic member component
// instead of specific to PIs?


export default function ProjectPI({
    data,
    showAvatar = true,
    showNameAndTitle = true,
    linkToProfile = false,
    variant = 'stroked',
}: {
    data: PrincipalInvestigator,
    showAvatar?: boolean,
    showNameAndTitle?: boolean,
    linkToProfile?: boolean
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

    const mainEl = (
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
                        height: '100%',
                    },
                    '& img': {
                        position: 'absolute',
                        top: 0,
                        left: 0,
                        opacity: imageLoaded ? 1 : 0,
                        transition: theme.transitions.create(
                            ['opacity'],
                            {
                                easing: theme.transitions.easing.swell,
                                duration: theme.transitions.duration.swell,
                            },
                        ),
                        width: '100%',
                        height: '100%',
                        objectFit: 'cover',
                        objectPosition: 'center',
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
                        color: 'ground.grade10',
                    },
                }}
            >
                <Typography
                    variant='pBodyBoldWt'
                >
                    {data.name}
                </Typography>

                {data.title &&
                    <Typography
                        variant='p3XS'
                    >
                        {data.title}
                    </Typography>
                }
            </Box>}
        </Box>
    )

    return (
        linkToProfile && data.profileUrl ? (
            <Link
                href={data.profileUrl}
            >
                {mainEl}
            </Link>
        ) : (
            mainEl
        )
    )
}
