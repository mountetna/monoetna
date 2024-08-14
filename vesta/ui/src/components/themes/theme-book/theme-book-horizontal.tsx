import * as React from 'react'
import Image from "next/image"
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import Collapse from '@mui/material/Collapse'
import Fade from '@mui/material/Fade'
import ButtonBase from '@mui/material/ButtonBase';

import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg'
import { ProjectCount, ThemeBookProps } from './shared';
import { useTheme } from '@mui/material';
import Link from '@/components/link/link';


export default function ThemeBookHorizontal({
    data,
    onClickSeeProjects,
    open,
    onSetOpen,
    onFinishOpen,
}: ThemeBookProps) {
    const theme = useTheme()

    const headingTextColor = data.textColor === 'light' ?
        theme.palette.utilityHighlight.main : theme.palette.ground.grade10

    const animationProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Box
            role='listitem'
            sx={() => {
                const styles = {
                    display: 'flex',
                    flexDirection: 'row',
                    // TODO: figure out why this is neccessary for proper animation
                    transform: 'translateY(0px)',
                    transition: theme.transitions.create(
                        ['transform'],
                        {
                            duration: theme.transitions.duration.quint,
                            easing: theme.transitions.easing.quint,
                        }
                    ),
                }

                return open ? styles : {
                    ...styles,
                    '&:hover, &:focus': {
                        transform: 'translateY(-30px)',
                    },
                }
            }}
        >
            {/* MAIN CONTENT */}
            <Collapse
                {...animationProps}
                onEntered={() => onFinishOpen()}
                orientation='horizontal'
            >
                <Fade
                    {...animationProps}
                    style={{
                        display: 'flex',
                    }}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            height: '100%',
                            p: '16px',
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'column',
                                flexGrow: '1',
                                justifyContent: 'space-between',
                                gap: '24px',
                            }}
                        >
                            <Box
                                sx={{
                                    display: 'flex',
                                    flexDirection: 'column',
                                    gap: '16px',
                                    textAlign: 'left',
                                    minWidth: '280px',
                                    [theme.breakpoints.up('desktop')]: {
                                        minWidth: '362px',
                                    },
                                }}
                            >
                                <Typography variant='h5'>
                                    {data.name}
                                </Typography>
                                <Typography variant='pLarge' component='div'>
                                    {data.description}
                                </Typography>
                            </Box>

                            {/* SEE PROJECTS BUTTON */}
                            <Box
                                sx={{
                                    display: 'flex',
                                    justifyContent: 'flex-end',
                                }}
                            >
                                <Link
                                    href={data.projectsLink}
                                    tabIndex={0}
                                    onClick={e => onClickSeeProjects(e, data.projectsLink)}
                                    sx={{
                                        display: 'inline-flex',
                                        alignItems: 'center',
                                        borderRadius: '60px',
                                        color: 'utilityHighlight.main',
                                        bgcolor: 'ground.grade10',
                                    }}
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
                                </Link>
                            </Box>
                        </Box>
                    </Box>
                </Fade>
            </Collapse>

            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'row',
                }}
            >
                {/* HEADING */}
                <ButtonBase
                    onClick={() => onSetOpen(!open)}
                    tabIndex={0}
                    sx={(theme) => ({
                        display: 'flex',
                        borderRadius: '20px',
                        borderTopRightRadius: open ? '0px' : '20px',
                        borderBottomRightRadius: open ? '0px' : '20px',
                        p: '16px',
                        bgcolor: open ? data.color : 'utilityWhite.main',
                        transition: theme.transitions.create(
                            ['background-color', 'border-radius'],
                            {
                                duration: theme.transitions.duration.quint,
                                easing: theme.transitions.easing.quint,
                            }
                        ),
                    })}
                >
                    <Box
                        sx={{
                            flexGrow: '1',
                            display: 'flex',
                            justifyContent: 'space-between',
                            flexDirection: 'column-reverse',
                            alignItems: 'center',
                            height: '100%',
                            rotate: '180deg',
                            '& .short-code, & .name': {
                                writingMode: 'vertical-lr',
                                textOrientation: 'mixed',
                            },
                        }}
                    >
                        <Typography
                            variant='h5'
                            component='span'
                            className='name'
                            sx={{
                                color: open ? headingTextColor : 'unset',
                                transition: theme.transitions.create(
                                    ['color'],
                                    {
                                        duration: theme.transitions.duration.quint,
                                        easing: theme.transitions.easing.quint,
                                    }
                                ),
                            }}
                        >
                            {data.name}
                        </Typography>
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'column',
                                alignItems: 'center',
                                gap: '10px',
                            }}
                        >
                            <Typography
                                variant='pSubtitleMonoCaps'
                                component='span'
                                className='short-code'
                                sx={{
                                    color: open ? headingTextColor : 'unset',
                                    transition: theme.transitions.create(
                                        ['color'],
                                        {
                                            duration: theme.transitions.duration.quint,
                                            easing: theme.transitions.easing.quint,
                                        }
                                    ),
                                }}
                            >
                                {data.name.slice(0, 4)}
                            </Typography>
                            <Typography
                                variant='pSubtitleMonoCaps'
                                component='div'
                                sx={{
                                    display: 'inline-flex',
                                }}
                            >
                                <ProjectCount
                                    count={data.projectCount}

                                />
                            </Typography>
                        </Box>
                    </Box>
                </ButtonBase>

                {/* COVER */}
                <Collapse
                    {...animationProps}
                    orientation='horizontal'

                >
                    <Fade
                        {...animationProps}
                        style={{
                            display: 'flex',
                            borderRadius: '20px',
                            borderTopLeftRadius: '0px',
                            borderBottomLeftRadius: '0px',
                            backgroundColor: data.color,
                        }}
                    >
                        <Box>
                            <Box
                                sx={(theme) => ({
                                    display: 'flex',
                                    overflow: 'hidden',
                                    p: '16px',
                                    width: '399px',
                                    height: '570px',
                                    [theme.breakpoints.up('desktop')]: {
                                        width: '516px',
                                        height: '671px',
                                    },
                                })}
                            >
                                <Image
                                    src={data.coverImage}
                                    alt={`${data.name} theme abstract image`}
                                    style={{
                                        width: '100%',
                                        height: '100%',
                                        objectFit: 'cover',
                                        objectPosition: 'center center',
                                        borderRadius: '16px',
                                    }}
                                />
                            </Box>
                        </Box>
                    </Fade>
                </Collapse>
            </Box>
        </Box >
    )
}