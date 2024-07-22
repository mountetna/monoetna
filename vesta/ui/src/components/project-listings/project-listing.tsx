'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography, { TypographyOwnProps } from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/material';
import ButtonBase from '@mui/material/ButtonBase';
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Image from 'next/image';
import { useSpring, animated, useIsomorphicLayoutEffect } from '@react-spring/web';
import { Group } from '@visx/group';
import { Pack, hierarchy } from '@visx/hierarchy';
import { useParentSize } from '@visx/responsive'
import _ from 'lodash';

import { ExternalProjectStatus, getExternalProjectStatus, Project } from "./models";
import { useWindowDimensions } from '@/lib/utils/responsive';
import ProjectStatus from './project-status';
import ProjectPI from './project-pi';

import projectListingBackshape from '/public/images/project-listing-backshape.svg'


interface ProjectAspect {
    title: string
    propName: keyof Project
    itemType: 'basic' | 'pill' | 'pi'
}

const projectAspects: ProjectAspect[] = [
    { title: 'Project Type', propName: 'type', itemType: 'basic', },
    { title: 'Data Types', propName: 'dataTypes', itemType: 'pill', },
    { title: 'Theme', propName: 'theme', itemType: 'pill', },
    { title: 'Principal Investigators', propName: 'principalInvestigators', itemType: 'pi', },
]

interface UserCountItem {
    radius: number
}


export default function ProjectListing({
    data,
    open,
    onSetOpen,
    onFinishOpen,
}: {
    data: Project,
    open: boolean,
    onSetOpen: (newState: boolean) => void,
    onFinishOpen: () => void,
}) {
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    // TODO
    const userHasAccess = false

    const { isResizing: isWindowResizing } = useWindowDimensions()

    // Manage lazy loading PI images
    const [isMainContentVisible, setIsMainContentVisible] = React.useState(false)
    React.useEffect(() => {
        if (open) {
            setIsMainContentVisible(true)
        }
    }, [open])

    // Manage open/close
    const mainContentRef = React.useRef<HTMLElement>()
    const [mainContentStyle, animateMainContentApi] = useSpring(() => ({
        height: '0px',
        opacity: 0,
        onRest: () => {
            if (open) {
                onFinishOpen()
            } else {
                setIsMainContentVisible(false)
            }
        },
    }), [open])

    useIsomorphicLayoutEffect(() => {
        animateMainContentApi.start(
            {
                height: `${open ? mainContentRef.current?.offsetHeight : 0}px`,
                opacity: open ? 1 : 0,
                config: {
                    easing: theme.transitions.easing.quintFn,
                    duration: theme.transitions.duration.quint,
                },
            }
        )
    }, [open, isWindowResizing])


    // Manage user count pack
    let {
        parentRef: userCountPackContainerRef,
        width: userCountPackContainerWidth,
        height: userCountPackContainerHeight,
    } = useParentSize({ debounceTime: 500, })
    // Add some variability to similar-looking packs
    const [packMarginPx] = React.useState<number>(_.random(3, 15, false))
    const [rotationDeg] = React.useState<number>(_.random(0, 360, false))

    userCountPackContainerWidth -= packMarginPx * 2
    userCountPackContainerHeight -= packMarginPx * 2

    // Make pack item size dependent on container width
    // so we can consistently size the center gap as a portion of container width
    const packItemAreaPx = Math.PI * ((userCountPackContainerWidth / 2) ** 2) / data.userCount
    const packItemRadiusPx = Math.sqrt(packItemAreaPx) / Math.PI
    const gapItemRadius = 0.33 * userCountPackContainerWidth / 2

    const packItems: UserCountItem[] = Array(data.userCount).fill({ radius: packItemRadiusPx })
    packItems.unshift({ radius: gapItemRadius })

    const pack = {
        children: packItems,
        name: 'root',
        radius: 0,
        distance: 0,
    }

    const packRoot = hierarchy<UserCountItem>(pack)
        .sum(d => d.radius ** 2)


    // Manage project aspects
    const projectAspectContentProps: TypographyOwnProps = {
        variant: 'pMedium',
        sx: {
            display: 'inline-flex',
            gap: '0px',
            color: 'ground.grade10',
        },
    }
    const projectAspectContentPillProps: TypographyOwnProps = {
        variant: 'pMediumMono',
        sx: {
            display: 'inline-flex',
            gap: '0px',
            color: 'ground.grade10',
            px: '10px',
            py: '2px',
            border: `1px solid ${theme.palette.ground.grade10}`,
            borderRadius: '20px',
        },
    }

    const allExtStatuses = Object.values(ExternalProjectStatus)
    const extStatus = getExternalProjectStatus(data)

    return (
        <Box
            className={`project-listing-item${open ? ' open' : ''}`}
            sx={{
                '&.open, &:hover, &:focus': {
                    boxShadow: '0px 8px 24px 0px #00000026'
                },
                transition: theme.transitions.create(
                    ['box-shadow'],
                    {
                        easing: theme.transitions.easing.ease,
                        duration: theme.transitions.duration.ease,
                    },
                ),
                position: 'relative',
                overflow: 'hidden',
            }}
        >
            <Box
                sx={{
                    position: 'absolute',
                    left: '0px',
                    top: '0px',
                    zIndex: 0,
                    opacity: isMainContentVisible ? 1 : 0,
                    [theme.breakpoints.up('tablet')]: {
                        left: '196px',
                        top: '15px',
                    },
                    transition: theme.transitions.create(
                        ['opacity'],
                        {
                            easing: theme.transitions.easing.ease,
                            duration: theme.transitions.duration.ease,
                        },
                    ),
                }}
            >
                <Image
                    src={projectListingBackshape}
                    alt='Abstract plant'
                />
            </Box>

            {/* HEADING */}
            <ButtonBase
                onClick={() => onSetOpen(!open)}
                tabIndex={0}
                sx={{
                    display: 'block',
                    textAlign: 'left',
                    width: '100%',
                    zIndex: 1,
                    position: 'relative',
                }}
            >
                <Box
                    sx={{
                        flexGrow: '1',
                        display: 'flex',
                        flexDirection: 'row',
                        alignItems: 'center',
                        justifyContent: 'space-between',
                        gap: '10px',
                        p: '16px',
                        // pb: open ? '12px' : '16px',
                        // transition: theme.transitions.create(
                        //     ['padding-bottom'],
                        //     {
                        //         easing: theme.transitions.easing.quint,
                        //         duration: theme.transitions.duration.quint,
                        //     },
                        // )
                    }}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            alignItems: 'center',
                            gap: '6px',
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                p: '6px',
                                bgcolor: data.theme.color,
                                borderRadius: '50%',
                            }}
                        >
                            <Image
                                src={data.theme.icon}
                                alt={`Abstract icon for "${data.theme.name}" theme`}
                                width={20}
                                height={20}
                                style={{
                                    backgroundColor: data.theme.color,
                                }}
                            />
                        </Box>

                        <Typography
                            variant={isMobile ? 'h6' : open ? 'h3' : 'h6'}
                            component='span'
                            sx={{
                                overflow: 'hidden',
                                textOverflow: 'ellipsis',
                                display: '-webkit-box',
                                WebkitLineClamp: 2,
                                WebkitBoxOrient: 'vertical',
                                opacity: isMobile && open ? 0 : 1,
                                transition: theme.transitions.create(
                                    ['opacity', 'font-weight', 'font-size'],
                                    {
                                        easing: theme.transitions.easing.quint,
                                        duration: theme.transitions.duration.quint,
                                    },
                                )
                            }}
                        >
                            {data.fullName}
                        </Typography>
                    </Box>

                    <Typography
                        variant='h6'
                        sx={{
                            borderRadius: '40px',
                            bgcolor: data.theme.color,
                            color: data.theme.textColor === 'light' ? 'utilityHighlight.main' : 'ground.grade10',
                            px: '14px',
                            py: '3px',
                            overflow: 'hidden',
                            textOverflow: 'ellipsis',
                            textAlign: 'center',
                            textWrap: 'nowrap',
                        }}
                    >
                        {data.theme.name}
                    </Typography>
                </Box>
            </ButtonBase>

            {/* MAIN CONTENT */}
            <animated.div
                style={{
                    overflow: 'hidden',
                    ...mainContentStyle,
                    position: 'relative',
                    zIndex: 1,
                }}
            >
                <Box
                    ref={mainContentRef}
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '18px',
                        p: '16px',
                        pt: '0px',
                        zIndex: 1,
                    }}
                >
                    <Typography
                        variant='h3MediumWt'
                        sx={{
                            [theme.breakpoints.up('tablet')]: {
                                display: 'none',
                            }
                        }}
                    >
                        {data.fullName}
                    </Typography>

                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '8px',
                            [theme.breakpoints.up('tablet')]: {
                                flexDirection: 'row',
                            },
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'column',
                                gap: '30px',
                                px: '8.21px',
                                [theme.breakpoints.up('tablet')]: {
                                    display: 'grid',
                                    gridTemplateColumns: 'repeat(6, 1fr)',
                                },
                                [theme.breakpoints.up('desktop')]: {
                                    gridTemplateColumns: 'repeat(11, 1fr)',
                                },
                            }}
                        >
                            {/* Image + Get Access */}
                            <Box
                                sx={{
                                    display: 'flex',
                                    flexDirection: 'column',
                                    gap: '24px',
                                    [theme.breakpoints.up('tablet')]: {
                                        gridColumn: 'span 3',
                                    },
                                    [theme.breakpoints.up('desktop')]: {
                                        gridColumn: 'span 5',
                                    },
                                }}
                            >
                                <Box
                                    sx={{
                                        position: 'relative',
                                        display: 'flex',
                                        aspectRatio: '344.57 / 247',
                                        borderRadius: '11.29px',
                                        overflow: 'hidden',
                                    }}
                                >
                                    {isMainContentVisible && <Image
                                        src={data.theme.imageComponents.projectBackground}
                                        alt={`Abstract background for ${data.theme.name} theme`}
                                        style={{
                                            width: '100%',
                                            height: '100%',
                                            objectFit: 'cover',
                                        }}
                                    />}
                                    <Box
                                        sx={{
                                            position: 'absolute',
                                            width: '100%',
                                            height: '100%',
                                            top: 0,
                                            left: 0,
                                        }}
                                    >
                                        <Typography
                                            variant='h5'
                                            component='div'
                                            sx={{
                                                width: '100%',
                                                height: '100%',
                                                display: 'flex',
                                                alignItems: 'center',
                                                justifyContent: 'center',
                                            }}
                                        >

                                            {/* Contributors count */}
                                            <Box
                                                sx={{
                                                    position: 'absolute',
                                                    left: 0,
                                                    top: 0,
                                                    width: '100%',
                                                    height: '100%',
                                                    display: 'flex',
                                                    alignItems: 'flex-end',
                                                    justifyContent: 'flex-end',
                                                    color: data.theme.textColor === 'light' ? 'utilityHighlight.main' : 'ground.grade10',
                                                    // color: data.theme.altColor,
                                                    px: '8px',
                                                    py: '5px',
                                                }}
                                            >
                                                <Box>
                                                    {`${data.userCount} Contributor${data.userCount === 1 ? '' : 's'}`}
                                                </Box>
                                            </Box>

                                            {/* User circles */}
                                            <Box
                                                sx={{
                                                    aspectRatio: '409 / 404',
                                                    height: '85%',
                                                    maxHeight: 'calc(100% - (1.5em * 2))',
                                                    borderRadius: '50%',
                                                    overflow: 'hidden',
                                                    position: 'relative',
                                                }}
                                            >
                                                {isMainContentVisible && <Image
                                                    src={data.theme.imageComponents.filtered}
                                                    alt={`Abstract foreground for ${data.theme.name} theme`}
                                                    style={{
                                                        width: '100%',
                                                        height: '100%',
                                                        objectFit: 'cover',
                                                    }}
                                                />}
                                                <Box
                                                    ref={userCountPackContainerRef}
                                                    sx={{
                                                        position: 'absolute',
                                                        width: '100%',
                                                        height: '100%',
                                                        left: 0,
                                                        top: 0,
                                                        display: 'flex',
                                                        justifyContent: 'center',
                                                        alignItems: 'center',
                                                    }}
                                                >
                                                    {isMainContentVisible && userCountPackContainerWidth > 10 &&
                                                        <svg
                                                            width={userCountPackContainerWidth}
                                                            height={userCountPackContainerHeight}
                                                            rx={userCountPackContainerWidth}
                                                            ry={userCountPackContainerHeight}
                                                            style={{
                                                                rotate: `${rotationDeg}deg`,
                                                            }}
                                                        >
                                                            <Pack<UserCountItem>
                                                                root={packRoot}
                                                                size={[userCountPackContainerWidth, userCountPackContainerHeight]}
                                                                padding={2}
                                                            >
                                                                {(packData) => {
                                                                    // skip root and center gap
                                                                    const circles = packData.descendants().slice(2)

                                                                    return (
                                                                        <Group>
                                                                            {circles.map((circle, i) => (
                                                                                <circle
                                                                                    key={`circle-${i}`}
                                                                                    r={circle.r}
                                                                                    cx={circle.x}
                                                                                    cy={circle.y}
                                                                                    fill={data.theme.altColor}
                                                                                    // onPointerEnter={e => e.currentTarget.style.fill = data.theme.color}
                                                                                    // onPointerLeave={e => e.currentTarget.style.fill = data.theme.altColor}
                                                                                    style={{
                                                                                        transition: theme.transitions.create(
                                                                                            ['fill'],
                                                                                            {
                                                                                                easing: theme.transitions.easing.ease,
                                                                                                duration: theme.transitions.duration.ease,
                                                                                            }
                                                                                        ),
                                                                                    }}
                                                                                />
                                                                            ))}
                                                                        </Group>
                                                                    )
                                                                }}
                                                            </Pack>
                                                        </svg>
                                                    }
                                                </Box>
                                            </Box>
                                        </Typography>
                                    </Box>
                                </Box>

                                <MUILink
                                    // TODO: href
                                    href={'#'}
                                    tabIndex={0}
                                    component={Link}
                                    underline='none'
                                    sx={{
                                        display: 'flex',
                                        minWidth: '100%',
                                        minHeight: '65.71px',
                                        px: '16.43px',
                                        py: '8.21px',
                                        borderRadius: '8.21px',
                                        alignItems: 'center',
                                        justifyContent: 'center',
                                        bgcolor: 'ground.grade10',
                                        transition: theme.transitions.create(
                                            ['background-color'],
                                            {
                                                easing: theme.transitions.easing.ease,
                                                duration: theme.transitions.duration.ease,
                                            },
                                        ),
                                        '&:hover, &:focus': {
                                            bgcolor: 'ground.grade25',
                                        },
                                    }}
                                >
                                    <Typography
                                        variant='pSubtitle'
                                        color='#F5F5F5'
                                    >
                                        {userHasAccess ? 'Visit Project' : 'Get Access'}
                                    </Typography>
                                </MUILink>
                            </Box>

                            {/* Description + Properties + Mobile/Desktop Status */}
                            <Box
                                sx={{
                                    display: 'flex',
                                    flexDirection: 'column',
                                    gap: '28.75px',
                                    [theme.breakpoints.up('tablet')]: {
                                        gridColumn: 'span 3',
                                    },
                                    [theme.breakpoints.up('desktop')]: {
                                        gridColumn: 'span 6',
                                    },
                                }}
                            >
                                <Typography variant='pMedium'>
                                    {data.description}
                                </Typography>

                                <Box
                                    sx={{
                                        display: 'flex',
                                        flexDirection: 'column',
                                        gap: '12.32px',
                                        '& .pill': {
                                            border: `1px solid ${theme.palette.ground.grade10}`,

                                        },
                                    }}
                                >
                                    {projectAspects.map(aspect => {
                                        return (
                                            <Box
                                                key={aspect.propName}
                                                sx={{
                                                    '& .pill-container, & .pi-container': {
                                                        display: 'flex',
                                                        flexWrap: 'wrap',
                                                        columnGap: '12.32px',
                                                        rowGap: '8px',
                                                        p: '4.11px',
                                                    },
                                                    '& .pi-container > *': {
                                                        // maxWidth: 'calc(100% / 4)',
                                                    },
                                                }}
                                            >
                                                <Typography
                                                    variant='pMediumMediumWt'
                                                    component='div'
                                                    sx={{
                                                        color: 'ground.grade10',
                                                        opacity: 0.4,
                                                        p: '4.11px',
                                                        pl: '0px',
                                                    }}
                                                >
                                                    {aspect.title}
                                                </Typography>

                                                <Box className={`${aspect.itemType}-container`}>
                                                    {['basic', 'pill'].includes(aspect.itemType) ? (

                                                        Array.isArray(data[aspect.propName]) ?
                                                            // @ts-ignore
                                                            data[aspect.propName].map((d, i) => (
                                                                <Typography
                                                                    key={`${aspect.propName}_${i}`}
                                                                    {...(aspect.itemType === 'basic' ? projectAspectContentProps : projectAspectContentPillProps)}
                                                                >
                                                                    {d}
                                                                </Typography>

                                                            ))
                                                            :
                                                            <Typography
                                                                {...(aspect.itemType === 'basic' ? projectAspectContentProps : projectAspectContentPillProps)}
                                                            >
                                                                {aspect.propName === 'theme' ? data.theme.name : data[aspect.propName]?.toString()}
                                                            </Typography>

                                                    ) : aspect.propName === 'principalInvestigators' ? (

                                                        data[aspect.propName].map((piData, i) => (
                                                            <ProjectPI
                                                                key={`${aspect.propName}_${i}`}
                                                                data={piData}
                                                                showAvatar={isMainContentVisible}
                                                            />
                                                        ))

                                                    ) : null
                                                    }
                                                </Box>
                                            </Box>
                                        )
                                    })}
                                </Box>

                                {/* Mobile + Desktop Status */}
                                <Box
                                    sx={{
                                        display: 'block',
                                        [theme.breakpoints.up('tablet')]: {
                                            display: 'none',
                                        },
                                        [theme.breakpoints.up('desktop')]: {
                                            display: 'block',
                                        },
                                    }}
                                >
                                    <ProjectStatus
                                        currentStatus={extStatus}
                                        allStatuses={allExtStatuses}
                                    />
                                </Box>
                            </Box>

                            {/* Tablet Status */}
                            <Box
                                sx={{
                                    display: 'none',
                                    [theme.breakpoints.up('tablet')]: {
                                        display: 'block',
                                        gridColumn: 'span 6'
                                    },
                                    [theme.breakpoints.up('desktop')]: {
                                        display: 'none',
                                    },
                                }}
                            >
                                <ProjectStatus
                                    currentStatus={extStatus}
                                    allStatuses={allExtStatuses}
                                />
                            </Box>
                        </Box>

                        {/* Dewey Decimal */}
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                gap: '8.63px',
                                p: '5.4px',
                                [theme.breakpoints.up('tablet')]: {
                                    flexDirection: 'column',
                                    alignItems: 'center',
                                    '& img': {
                                        rotate: '90deg',
                                    },
                                },
                            }}
                        >
                            <Image
                                src={data.theme.icon}
                                alt={`Abstract icon for "${data.theme.name}" theme`}
                                width={22.6}
                                height={22.6}
                            />

                            <Typography
                                variant='pBodyMono'
                                sx={{
                                    display: 'flex',
                                    flexDirection: 'row',
                                    gap: '3.4px',
                                    textTransform: 'uppercase',
                                    [theme.breakpoints.up('tablet')]: {
                                        writingMode: 'vertical-lr',
                                        textOrientation: 'mixed',
                                    },
                                }}
                            >
                                <Box>
                                    {`${data.theme.name.slice(0, 2)}.${data.name}`}
                                </Box>
                                <Box>
                                    {`${data.type.slice(0, 4)}.${data.status.slice(0, 4)}`}
                                </Box>
                                <Box>
                                    {data.startDate.getFullYear()}
                                </Box>
                            </Typography>
                        </Box>
                    </Box>
                </Box>
            </animated.div>
        </Box>
    )
}