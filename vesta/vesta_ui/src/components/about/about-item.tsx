'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import { containerPadding } from '@/theme'
import Image, { StaticImageData } from 'next/image'


export interface Link {
    header?: string
    blurb?: string
    href: string
    label: string
}

export interface ImageProps {
    src: StaticImageData
    alt: string
}

export default function AboutItem({
    header,
    body,
    link,
    image,
}: {
    header: string,
    body: string,
    link?: Link,
    image: ImageProps,
}) {
    return (
        <Box
            className='about-item'
            sx={{
                width: '100vw',
                maxWidth: 'desktopLg',
            }}
        >
            <Box
                sx={{
                    ...containerPadding,
                }}
            >
                <Box
                    sx={(theme) => ({
                        bgcolor: 'utilityWhite.main',
                        borderRadius: '30px',
                        px: '8px',
                        py: '8px',
                        [theme.breakpoints.up('tablet')]: {
                            px: '30px',
                            py: '30px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            px: '30px',
                            py: '30px',
                            display: 'flex',
                            flexDirection: 'row-reverse',
                            columnGap: '34px',
                        },
                    })}
                >
                    <Box
                        sx={(theme) => ({
                            aspectRatio: '649 / 601',
                            mb: '16px',
                            [theme.breakpoints.up('desktop')]: {
                                mb: '0',
                                width: '45%',
                            },
                            '& img': {
                                borderRadius: '30px',
                                objectFit: 'cover',
                                width: '100%',
                                height: '100%',
                            },
                        })}
                    >
                        <Image
                            src={image.src}
                            alt={image.alt}
                        />
                    </Box>
                    <Box
                        sx={(theme) => ({
                            [theme.breakpoints.up('desktop')]: {
                                width: '55%',
                                display: 'flex',
                                flexDirection: 'column',
                                justifyContent: 'space-between',
                            },
                        })}
                    >
                        <Box
                            sx={(theme) => ({
                                mb: '32px',
                            })}
                        >
                            <Typography
                                variant='h2'
                                sx={(theme) => ({
                                    mb: '24px',
                                })}
                            >
                                {header}
                            </Typography>
                            <Typography variant='pLarge'>
                                {body}
                            </Typography>
                        </Box>
                        {link && <Box>
                            {(link.header || link.blurb) &&
                                <Box
                                    sx={{
                                        mb: '24px',
                                    }}
                                >
                                    {link.header && <Typography
                                        variant='h6BoldWt'
                                        component='h6'
                                    >
                                        {link.header}
                                    </Typography>}
                                    {link.blurb && <Typography
                                        variant='pLarge'
                                        sx={(theme) => ({
                                            mb: '24px',
                                        })}
                                    >
                                        {link.blurb}
                                    </Typography>}
                                </Box>}
                            <Box
                                sx={(theme) => ({
                                    p: '8px',
                                    [theme.breakpoints.up('tablet')]: {
                                        p: '0',
                                    },
                                })}
                            >
                                <MUILink
                                    href={link.href}
                                    tabIndex={0}
                                    component={Link}
                                    underline='none'
                                >
                                    <Typography
                                        variant='h5'
                                        sx={{
                                            display: 'inline-block',
                                            px: '16px',
                                            py: '8px',
                                            borderRadius: '8px',
                                            bgcolor: 'ground.grade10',
                                            color: '#F5F5F5',
                                        }}
                                    >
                                        {link.label}
                                    </Typography>
                                </MUILink>
                            </Box>
                        </Box>}
                    </Box>
                </Box>
            </Box>
        </Box>
    )
}