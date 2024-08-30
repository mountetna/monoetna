'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import ButtonBase from '@mui/material/ButtonBase';
import Collapse from '@mui/material/Collapse'
import Fade from '@mui/material/Fade'
import { TransitionProps } from '@mui/material/transitions';
import _ from 'lodash';

import { ExternalProjectStatus, getExternalProjectStatus, Project, ProjectHeadingInfoSet } from './models';
import ProjectStatus from './project-status';
import ProjectPI from './project-pi';
import projectListingBackshape from '/public/images/project-listing-backshape.svg'
import ProjectHeadingInfo from './project-heading-info';
import Image from '../image/image';
import { useBreakpoint } from '@/lib/utils/responsive';
import Pill from '../pill/pill';
import ProjectUserCountPack from './project-user-count-pack';
import Link from '../link/link';


export enum Classes {
    iconAndHeadingContainer = 'icon-heading-container'
}

interface ProjectAspectBase {
    title: string
    kind: string
}

interface ProjectAspect extends ProjectAspectBase {
    title: string
    propName: keyof Project
    itemType: 'basic' | 'pill' | 'pi' | 'counted'
    kind: 'default'
}

interface CountedProp {
    label: string
    propName: keyof Project
}

interface CountedProjectAspect extends ProjectAspectBase {
    title: string
    props: CountedProp[]
    kind: 'counted'
}

const projectAspects: (ProjectAspect | CountedProjectAspect)[] = [
    { title: 'Project Type', propName: 'type', itemType: 'basic', kind: 'default' },
    { title: 'Lab Metrics', props: [{ label: 'Samples', propName: 'sampleCount' }, { label: 'Assays', propName: 'assayCount' }], kind: 'counted' },
    { title: 'Data Types', propName: 'dataTypes', itemType: 'pill', kind: 'default' },
    { title: 'Clinical Data', propName: 'hasClinicalData', itemType: 'basic', kind: 'default' },
    { title: 'Theme', propName: 'theme', itemType: 'pill', kind: 'default' },
    { title: 'Principal Investigators', propName: 'principalInvestigators', itemType: 'pi', kind: 'default' },
]


export default function ProjectListing({
    data,
    open,
    headingInfoSet,
    onSetOpen,
    onFinishOpen,
}: {
    data: Project,
    open: boolean,
    headingInfoSet: ProjectHeadingInfoSet,
    onSetOpen: (newState: boolean) => void,
    onFinishOpen: () => void,
}) {
    const theme = useTheme()
    const breakpoint = useBreakpoint()
    const isMobile = breakpoint === 'mobile'
    const isDesktop = (['desktop', 'desktopLg'] as (typeof breakpoint)[]).includes(breakpoint)

    // TODO
    const userHasAccess = false

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    // Manage lazy loading
    const [hasOpened, setHasOpened] = React.useState(false)
    React.useEffect(() => {
        if (open) {
            setHasOpened(true)
        }
    }, [open])

    const allExtStatuses = Object.values(ExternalProjectStatus)
    const extStatus = getExternalProjectStatus(data)

    return (
        <Box
            className={`project-listing-item${open ? ' open' : ''}`}
            sx={{
                '&.open, &:hover, &:focus': {
                    boxShadow: '0px 8px 24px 0px #00000026',
                },
                '&.open': {

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
                    left: '216px',
                    top: '35px',
                    zIndex: 0,
                    opacity: open ? 1 : 0,
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
                    }}
                >
                    <Box
                        className={Classes.iconAndHeadingContainer}
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            alignItems: 'center',
                            gap: '6px',
                            [theme.breakpoints.up('tablet')]: {
                                gap: '8px'
                            },
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
                                hideBeforeLoad={false}
                            />
                        </Box>

                        <Typography
                            variant={isMobile ? 'h6' : open ? 'h3' : 'h6'}
                            component='span'
                            sx={{
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

                    <ProjectHeadingInfo
                        projectData={data}
                        projectOpen={open}
                        infoSet={headingInfoSet}
                        variant={isDesktop ? 'default' : 'small'}
                    />
                </Box>
            </ButtonBase>

            {/* MAIN CONTENT */}
            <Collapse
                {...animationProps}
                onEntered={() => onFinishOpen()}
                sx={{
                    position: 'relative',
                    zIndex: 1,
                }}
            >
                <Fade
                    {...animationProps}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '18px',
                            p: '16px',
                            pt: '0px',
                            zIndex: 1,
                            [theme.breakpoints.up('tablet')]: {
                                gap: '0px',
                                pt: '36px',
                            },
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
                                pt: '20px',
                                [theme.breakpoints.up('tablet')]: {
                                    flexDirection: 'row',
                                    gap: '16px',
                                    pb: '60px',
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
                                        columnGap: '30px',
                                        rowGap: '16px',
                                    },
                                    [theme.breakpoints.up('desktop')]: {
                                        gridTemplateColumns: 'repeat(11, 1fr)',
                                        columnGap: '32px',
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
                                        {hasOpened && <Image
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
                                                        px: '8px',
                                                        py: '5px',
                                                    }}
                                                >
                                                    <Box>
                                                        {`${data.userCount} User${data.userCount === 1 ? '' : 's'}`}
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
                                                    {hasOpened && (
                                                        <Image
                                                            src={data.theme.imageComponents.filtered}
                                                            alt={`Abstract foreground for ${data.theme.name} theme`}
                                                            style={{
                                                                width: '100%',
                                                                height: '100%',
                                                                objectFit: 'cover',
                                                            }}
                                                        />
                                                    )}

                                                    {hasOpened && (
                                                        <ProjectUserCountPack
                                                            userCount={data.userCount}
                                                            color={data.theme.altColor}
                                                            visible={open}
                                                        />
                                                    )}
                                                </Box>
                                            </Typography>
                                        </Box>
                                    </Box>

                                    <Link
                                        // TODO: href
                                        href={'#'}
                                        tabIndex={0}
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
                                    </Link>
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
                                                    key={aspect.title}
                                                    sx={{
                                                        [theme.breakpoints.up('desktop')]: {
                                                            display: 'grid',
                                                            gridTemplateColumns: 'repeat(3, 1fr)',
                                                            columnGap: '53.39px',
                                                        },
                                                        '& .pill-container, & .pi-container': {
                                                            display: 'flex',
                                                            flexWrap: 'wrap',
                                                            columnGap: '12.32px',
                                                            rowGap: '8px',
                                                            p: '4.11px',
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
                                                            [theme.breakpoints.up('desktop')]: {
                                                                gridColumn: 'span 1',
                                                            },
                                                        }}
                                                    >
                                                        {aspect.title}
                                                    </Typography>

                                                    {aspect.kind === 'default' ? (
                                                        <Box
                                                            className={`${aspect.itemType}-container`}
                                                            sx={{
                                                                [theme.breakpoints.up('desktop')]: {
                                                                    gridColumn: 'span 2',
                                                                },
                                                            }}
                                                        >
                                                            {['basic', 'pill'].includes(aspect.itemType) ? (

                                                                Array.isArray(data[aspect.propName]) ?
                                                                    // @ts-ignore
                                                                    data[aspect.propName].map((d, i) => {
                                                                        const key = `${aspect.propName}_${i}`

                                                                        if (aspect.itemType === 'basic') {
                                                                            return (
                                                                                <Typography
                                                                                    key={key}
                                                                                    variant='pMedium'
                                                                                    sx={{
                                                                                        gap: '0px',
                                                                                        color: 'ground.grade10',
                                                                                    }}
                                                                                >
                                                                                    {d}
                                                                                </Typography>
                                                                            )
                                                                        } else {
                                                                            return (
                                                                                <Pill
                                                                                    key={key}
                                                                                    // @ts-ignore
                                                                                    label={d}
                                                                                    typographyVariant='pMediumMono'
                                                                                    variant='stroked'
                                                                                    sx={{
                                                                                        textAlign: 'left'
                                                                                    }}
                                                                                />
                                                                            )
                                                                        }

                                                                    })
                                                                    : (
                                                                        aspect.itemType === 'basic' ? (
                                                                            <Typography
                                                                                key={aspect.propName}
                                                                                variant='pMedium'
                                                                                sx={{
                                                                                    gap: '0px',
                                                                                    color: 'ground.grade10',
                                                                                }}
                                                                            >
                                                                                {aspect.propName === 'theme' ? data.theme.name : data[aspect.propName]?.toString()}
                                                                            </Typography>
                                                                        ) : (
                                                                            <Pill
                                                                                key={aspect.propName}
                                                                                // @ts-ignore
                                                                                label={aspect.propName === 'theme' ? data.theme.name : data[aspect.propName]?.toString()}
                                                                                typographyVariant='pMediumMono'
                                                                                variant='stroked'
                                                                                sx={{
                                                                                    textAlign: 'left'
                                                                                }}
                                                                            />
                                                                        )
                                                                    )


                                                            ) : aspect.propName === 'principalInvestigators' ? (

                                                                data[aspect.propName].map((piData, i) => (
                                                                    <ProjectPI
                                                                        key={`${aspect.propName}_${i}`}
                                                                        data={piData}
                                                                        showAvatar={hasOpened}
                                                                    />
                                                                ))

                                                            ) : null
                                                            }
                                                        </Box>
                                                    ) : (aspect.kind === 'counted' ? (
                                                        <Box
                                                            className='pill-container'
                                                            sx={{
                                                                [theme.breakpoints.up('desktop')]: {
                                                                    gridColumn: 'span 2',
                                                                },
                                                            }}
                                                        >
                                                            {aspect.props.map(prop => (
                                                                <Pill
                                                                    key={prop.label}
                                                                    label={prop.label}
                                                                    // @ts-ignore
                                                                    count={data[prop.propName]}
                                                                    typographyVariant='pMediumMono'
                                                                    variant='counter'
                                                                    sx={{
                                                                        textAlign: 'left'
                                                                    }}
                                                                />
                                                            ))}
                                                        </Box>
                                                    ) : null
                                                    )}
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
                                            variant={isMobile ? 'compact' : 'default'}
                                            visible={open}
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
                                        variant='default'
                                        visible={open}
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
                                    },
                                    '& img': {
                                        [theme.breakpoints.up('tablet')]: {
                                            rotate: '90deg',
                                            position: 'relative',
                                            left: '0.6px',
                                        },
                                        [theme.breakpoints.up('desktop')]: {
                                            left: '0.1px',
                                        },
                                    },
                                }}
                            >
                                <Image
                                    src={data.theme.icon}
                                    alt={`Abstract icon for "${data.theme.name}" theme`}
                                    width={isDesktop ? 28 : 22.6}
                                    height={isDesktop ? 28 : 22.6}
                                />

                                <Typography
                                    variant={isDesktop ? 'pTitleMono' : 'pBodyMono'}
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
                                        {`${data.theme.name.slice(0, 4)}.${data.name}`}
                                    </Box>
                                    <Box>
                                        {`${data.type.slice(0, 4)}.${getExternalProjectStatus(data).slice(0, 4)}`}
                                    </Box>
                                    <Box>
                                        {data.startDate.getFullYear()}
                                    </Box>
                                </Typography>
                            </Box>
                        </Box>
                    </Box>
                </Fade>
            </Collapse>
        </Box>
    )
}