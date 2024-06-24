import * as React from 'react'
import Image from "next/image"
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import { useSpring, animated } from 'react-spring';
import ButtonBase from '@mui/material/ButtonBase';

import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg'
import { ProjectCount, ThemeData } from './shared';


export default function ThemeBookVertical({
    data,
    onClickSeeProjects,
    open,
    onSetOpen,
    onFinishOpen,
}: {
    data: ThemeData,
    onClickSeeProjects: (event: React.MouseEvent<HTMLAnchorElement>, href: string) => void,
    open: boolean,
    onSetOpen: (newState: boolean) => void,
    onFinishOpen: () => void,
}) {
    const mainContentRef = React.useRef<HTMLElement>()
    const [mainContentStyle, animateMainContentApi] = useSpring(() => {
        return {
            height: '0px',
            opacity: 0,
            // TODO: fix. why doesn't this trigger?
            onRest: () => onFinishOpen,
        }
    }, [open])

    React.useEffect(() => {
        animateMainContentApi.start({
            height: `${open ? mainContentRef.current?.offsetHeight : 0}px`,
            opacity: open ? 1 : 0,
        })
    }, [open])

    return (
        <Box
            role='listitem'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                p: '16px',
                bgcolor: 'utilityWhite.main',
                borderRadius: '20px',
            })}
        >
            {/* HEADING */}
            <ButtonBase
                onClick={() => onSetOpen(!open)}
                tabIndex={0}
            >
                <Box
                    sx={(theme) => ({
                        display: 'flex',
                        justifyContent: 'space-between',
                        flexGrow: '1',
                    })}
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
            <animated.div
                style={{
                    overflow: 'hidden',
                    ...mainContentStyle,
                }}
            >
                <Box
                    ref={mainContentRef}
                    sx={(theme) => ({
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '19px',
                        pt: '24px',
                    })}
                >
                    <Box
                        sx={(theme) => ({
                            display: 'flex',
                            justifyContent: 'space-between',
                            gap: '16px',
                            p: '16px',
                            borderRadius: '20px',
                            bgcolor: data.color,
                            objectFit: 'contain',
                        })}
                    >
                        <Typography
                            variant='pMediumMonoCaps'
                            component='div'
                            sx={(theme) => ({
                                alignSelf: 'flex-end',
                                display: 'flex',
                                [theme.breakpoints.up('tablet')]: {
                                    display: 'none',
                                },
                            })}
                        >
                            <ProjectCount
                                count={data.project_count}
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
                                src={data.image}
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
                                href={data.projects_link}
                                tabIndex={0}
                                component={Link}
                                underline='none'
                                color='utilityHighlight.main'
                                bgcolor='ground.grade10'
                                onClick={e => onClickSeeProjects(e, data.projects_link)}
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
            </animated.div>
        </Box>
    )
}