'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';

import { PrincipalInvestigator } from "./models";
import Image from 'next/image';


// Maybe make this a generic member component
// instead of specific to PIs?

// TODO: add profile link if available
export default function ProjectPI({
    data,
    showAvatar = true,
}: {
    data: PrincipalInvestigator,
    showAvatar?: boolean,
}) {
    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'row',
                gap: '5.13px',
            }}
        >
            <Box
                sx={{
                    width: '41px',
                    height: '41px',
                    borderRadius: '50%',
                    overflow: 'hidden',
                    '& > *': {
                        width: '100%',
                        minWidth: '100%',
                        height: '100%',
                        minHeight: '100%',
                    },
                }}
            >
                {showAvatar && data.imageUrl ?
                    <Image
                        src={data.imageUrl}
                        alt={`Profile image for ${data.name}`}
                        width={41}
                        height={41}
                    />
                    :
                    <Box
                        sx={{
                            display: 'flex',
                            justifyContent: 'center',
                            alignItems: 'center',
                            bgcolor: data.color,
                        }}
                    >
                        <Typography variant='pBodyBoldWt'>
                            {data.name[0].toUpperCase()}
                        </Typography>
                    </Box>
                }
            </Box>

            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    // justifyContent: 'space-around',
                    '& *': {
                        textOverflow: 'ellipsis',
                        textWrap: 'nowrap',
                        display: '-webkit-box',
                        WebkitLineClamp: 1,
                        WebkitBoxOrient: 'vertical',
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
            </Box>
        </Box>
    )
}