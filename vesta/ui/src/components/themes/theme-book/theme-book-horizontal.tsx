import * as React from 'react'
import Image from "next/image"
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import { useSpring, animated, easings } from '@react-spring/web';
import ButtonBase from '@mui/material/ButtonBase';

import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg'
import { ProjectCount, ThemeBookProps } from './shared';
import { useWindowDimensions } from '@/lib/utils/responsive';
import { useTheme } from '@mui/material';


export default function ThemeBookHorizontal({
    data,
    onClickSeeProjects,
    open,
    onSetOpen,
    onFinishOpen,
}: ThemeBookProps) {
    const theme = useTheme()

    const mainContentRef = React.useRef<HTMLElement>()
    const coverRef = React.useRef<HTMLElement>()

    const springParams = {
        width: '0px',
        opacity: 0,
        onRest: () => open && onFinishOpen(),
    }
    const [mainContentStyle, animateMainContentApi] = useSpring(() => ({ ...springParams }), [])
    const [coverStyle, animateCoverApi] = useSpring(() => ({ ...springParams }), [])

    const { isResizing: isWindowResizing } = useWindowDimensions()
    React.useEffect(() => {
        animateMainContentApi.start({
            width: `${open ? mainContentRef.current?.offsetWidth : 0}px`,
            opacity: open ? 1 : 0,
            config: {
                easing: theme.transitions.easing.quintFn,
                duration: theme.transitions.duration.quint,
            },
        })
        animateCoverApi.start({
            width: `${open ? coverRef.current?.offsetWidth : 0}px`,
            opacity: open ? 1 : 0,
            config: {
                easing: theme.transitions.easing.quintFn,
                duration: theme.transitions.duration.quint,
            },
        })
    }, [open, isWindowResizing])

    const headingTextColor = data.textColor === 'light' ?
        theme.palette.utilityHighlight.main : theme.palette.ground.grade10

    return (
        <Box
            role='listitem'
            sx={() => {
                const styles = {
                    display: 'flex',
                    flexDirection: 'row',
                    transition: theme.transitions.create(
                        ['transform'],
                        {
                            duration: theme.transitions.duration.ease,
                            easing: theme.transitions.easing.ease,
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
            <animated.div
                style={{
                    display: 'flex',
                    overflow: 'hidden',
                    ...mainContentStyle,
                }}
            >
                <Box
                    ref={mainContentRef}
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
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
                            <MUILink
                                href={data.projects_link}
                                tabIndex={0}
                                component={Link}
                                underline='none'
                                onClick={e => onClickSeeProjects(e, data.projects_link)}
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
                            </MUILink>
                        </Box>
                    </Box>
                </Box>
            </animated.div>

            {/* HEADING */}
            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'row',
                }}
            >
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
                                duration: theme.transitions.duration.ease,
                                easing: theme.transitions.easing.ease,
                            }
                        ),
                    })}
                >
                    <Box
                        sx={{
                            flexGrow: '1',
                            display: 'flex',
                            justifyContent: 'space-between',
                            flexDirection: 'column',
                            alignItems: 'center',
                            height: '100%',
                            '& .short-code, & .name': {
                                writingMode: 'vertical-lr',
                                textOrientation: 'mixed',
                                rotate: '180deg',
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
                                        duration: theme.transitions.duration.ease,
                                        easing: theme.transitions.easing.ease,
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
                                gap: '10px',
                            }}
                        >
                            <Typography
                                variant='pSubtitleMonoCaps'
                                component='div'
                                sx={{
                                    display: 'inline-flex',
                                }}
                            >
                                <ProjectCount
                                    count={data.project_count}

                                />
                            </Typography>
                            <Typography
                                variant='pSubtitleMonoCaps'
                                component='span'
                                className='short-code'
                                sx={{
                                    color: open ? headingTextColor : 'unset',
                                    transition: theme.transitions.create(
                                        ['color'],
                                        {
                                            duration: theme.transitions.duration.ease,
                                            easing: theme.transitions.easing.ease,
                                        }
                                    ),
                                }}
                            >
                                {data.name.slice(0, 4)}
                            </Typography>
                        </Box>
                    </Box>
                </ButtonBase>

                {/* COVER */}
                <animated.div
                    style={{
                        display: 'flex',
                        borderRadius: '20px',
                        borderTopLeftRadius: '0px',
                        borderBottomLeftRadius: '0px',
                        backgroundColor: data.color,
                        overflow: 'hidden',
                        ...coverStyle,
                    }}
                >
                    <Box
                        ref={coverRef}
                    >
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
                                src={data.image}
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
                </animated.div>
            </Box>
        </Box>
    )
}