import * as React from 'react'
import Image from "next/image"
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import MUILink from '@mui/material/Link';
import Collapse from '@mui/material/Collapse'
import Fade from '@mui/material/Fade'
import { TransitionProps } from '@mui/material/transitions';
import Link from 'next/link'
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material';

import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg'
import { ProjectCount, ThemeBookProps } from './shared';


export default function ThemeBookVertical({
    data,
    onClickSeeProjects,
    open,
    onSetOpen,
    onFinishOpen,
}: ThemeBookProps) {
    const theme = useTheme()

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Box
            role='listitem'
            sx={{
                display: 'flex',
                flexDirection: 'column',
                bgcolor: 'utilityWhite.main',
                borderRadius: '20px',
            }}
        >
            {/* HEADING */}
            <ButtonBase
                onClick={() => onSetOpen(!open)}
                tabIndex={0}
                sx={{
                    p: '16px',
                }}
            >
                <Box
                    sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        flexGrow: '1',
                    }}
                >
                    <Typography
                        variant='pSubtitleMonoCaps'
                        component='span'
                    >
                        {data.name.slice(0, 4)}
                    </Typography>
                    <Typography
                        variant='h5'
                        component='span'
                    >
                        {data.name}
                    </Typography>
                </Box>
            </ButtonBase>

            {/* MAIN CONTENT */}
            <Collapse
                {...animationProps}
                onEntered={() => onFinishOpen()}
            >
                <Fade
                    {...animationProps}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '19px',
                            p: '8px 16px 16px 16px',
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                justifyContent: 'space-between',
                                gap: '16px',
                                p: '16px',
                                borderRadius: '20px',
                                bgcolor: data.color,
                                objectFit: 'contain',
                            }}
                        >
                            <Typography
                                variant='pMediumMonoCaps'
                                component='div'
                                sx={(theme) => ({
                                    alignSelf: 'flex-end',
                                    display: 'flex',
                                    rotate: '180deg',
                                    [theme.breakpoints.up('tablet')]: {
                                        display: 'none',
                                    },
                                })}
                            >
                                <ProjectCount
                                    count={data.projectCount}
                                />
                            </Typography>
                            <Box
                                sx={(theme) => ({
                                    display: 'flex',
                                    alignSelf: 'flex-start',
                                    '& img': {
                                        width: '100%',
                                        minHeight: '272px',
                                        height: '100%',
                                        aspectRatio: '270 / 272',
                                        objectFit: 'cover',
                                        objectPosition: 'center center',
                                        borderRadius: '16px',
                                    },
                                })}
                            >
                                <Image
                                    src={data.coverImage}
                                    alt={`${data.name} theme abstract image`}
                                />
                            </Box>
                        </Box>

                        <Box
                            sx={(theme) => ({
                                display: 'flex',
                                flexDirection: 'column',
                                gap: '24px',
                            })}
                        >
                            <Box
                                sx={(theme) => ({
                                    display: 'flex',
                                    flexDirection: 'column',
                                    gap: '16px',
                                    textAlign: 'left',
                                })}
                            >
                                <Typography variant='h5'>
                                    {data.name}
                                </Typography>
                                <Typography variant='pLarge'>
                                    {data.description}
                                </Typography>
                            </Box>

                            <Box
                                sx={{
                                    display: 'flex',
                                    justifyContent: 'flex-end',
                                }}
                            >
                                <MUILink
                                    href={data.projectsLink}
                                    tabIndex={0}
                                    component={Link}
                                    underline='none'
                                    color='utilityHighlight.main'
                                    bgcolor='ground.grade10'
                                    onClick={e => onClickSeeProjects(e, data.projectsLink)}
                                    sx={(theme) => ({
                                        display: 'inline-flex',
                                        alignItems: 'center',
                                        borderRadius: '60px',
                                    })}
                                >
                                    <Typography
                                        variant='pLarge'
                                        sx={{
                                            display: 'inline-flex',
                                            p: '10px 10px 10px 24px',
                                        }}
                                    >
                                        See projects
                                    </Typography>
                                    <Box
                                        sx={{
                                            display: 'flex',
                                            p: '12px',
                                        }}
                                    >
                                        <Image
                                            src={arrowUpRightLight}
                                            alt='Arrow pointing up-right'
                                            width={40}
                                            height={40}
                                        />
                                    </Box>
                                </MUILink>
                            </Box>
                        </Box>
                    </Box>
                </Fade>
            </Collapse>
        </Box>
    )
}