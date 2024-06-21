import * as React from 'react'
import Image, { StaticImageData } from "next/image"
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import { useSpring, animated } from 'react-spring';
import ButtonBase from '@mui/material/ButtonBase';

import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg'
import { useMediaQuery, useTheme } from '@mui/material';


export interface ThemeData {
    name: string
    description: string
    project_count: number
    projects_link: string
    color: string
    image: StaticImageData
}


export default function ThemeBook({
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
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    const mainContentRef = React.useRef<HTMLElement>()
    const [mainContentStyle, animateMainContentApi] = useSpring(
        () => ({
            height: '0px',
            opacity: 0,
            // TODO: fix. why doesn't this trigger?
            // onRest: () => onFinishOpen(),
        }),
        []
    )

    React.useEffect(() => {
        animateMainContentApi.start({
            height: `${open ? mainContentRef.current?.offsetHeight : 0}px`,
            opacity: open ? 1 : 0,
        })
    }, [open])

    const content = (
        <Box
            role='listitem'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                // rowGap: '24px',
                p: '16px',
                bgcolor: 'utilityWhite.main',
                borderRadius: '20px',
            })}
        >
            {/* HEADING */}
            <Box
                sx={(theme) => ({
                    display: 'flex',
                    justifyContent: 'space-between',
                })}
            >
                <Typography
                    variant='pSubtitleMonoCaps'
                    component='span'
                >
                    {data.name.slice(0, 4)}
                </Typography>
                <Typography
                    variant='pSubtitleMonoCaps'
                    component='div'
                    sx={(theme) => ({
                        display: 'none',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'inline-block',
                        },
                    })}
                >
                    <ProjectCount
                        count={data.project_count}

                    />
                </Typography>
                <Typography
                    variant='h5'
                    component='span'
                >
                    {data.name}
                </Typography>
            </Box>

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

    return (
        isMobile ?
            <ButtonBase
                onClick={() => onSetOpen(!open)}
                tabIndex={0}
            >
                {content}
            </ButtonBase>
            : content
    )
}


function ProjectCount({
    count,
}: {
    count: number,
}) {
    return (
        <Box
            component='span'
            sx={{
                display: 'inline-flex',
                columnGap: '10px',
                py: '10px',
                borderRadius: '30px',
                color: 'utilityWhite.main',
                bgcolor: 'ground.grade10',
                writingMode: 'vertical-lr',
                textOrientation: 'mixed',
                rotate: '180deg',
            }}
        >
            <Box component='span'>
                {count}
            </Box>
            <Box component='span'>
                {`PROJECT${count > 1 ? 'S' : ''}`}
            </Box>
        </Box>
    )
}