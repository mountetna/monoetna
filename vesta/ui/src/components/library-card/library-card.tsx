import * as React from 'react';
import Box from '@mui/system/Box'
import { SxProps, Typography, useTheme } from '@mui/material';
import Image from 'next/image';

import { User } from '../user/models';

import logoDarkSrc from '/public/images/logo/logo-dark.svg'
import { useDeviceOrientation, usePointerPosition } from '@/lib/utils/position';


export enum Classes {
    root = 'library-card',
    side = 'library-card-side',
    depthLayer = 'library-card-depth-layer'
}


interface Stat {
    key: keyof User
    label: string
}

const STATS: Stat[] = [
    // { key: 'contributions', label: 'Contributions', },
    { key: 'projectMemberships', label: 'Projects', },
]

interface Props {
    user: User,
}


function LibraryCard(props: Props, ref: React.ForwardedRef<unknown>) {
    const theme = useTheme()
    const user = props.user

    const depth = 5
    const outlineSize = 0.5
    const rotationIntensity = 15
    const shineSize = 1500

    // Manage effects
    // - Orientation
    const pointerCoords = usePointerPosition('fraction')
    const deviceOrientationCoords = useDeviceOrientation()
    let shineX = 0
    let shineY = 0
    let transform = 'unset'

    if (deviceOrientationCoords !== null) {
        const xCoord = deviceOrientationCoords.x
        const yCoord = deviceOrientationCoords.y
        const zCoord = deviceOrientationCoords.z
        // const rotX = (yCoord - 0.5) / 0.5
        // const rotY = -(xCoord - 0.5) / 0.5
        // const rotMagnitude = Math.max(Math.abs(rotX), Math.abs(rotY))
        const rotDeg = Math.round(rotationIntensity)
        // const rotDeg = Math.round(rotationIntensity * rotMagnitude)
        transform = `rotateZ(${zCoord}deg) rotateX(${xCoord}deg) rotateY(${yCoord}deg)`

        shineX = xCoord
        shineY = yCoord
    } else if (pointerCoords !== null) {
        const xCoord = pointerCoords.x
        const yCoord = pointerCoords.y
        const rotX = (yCoord - 0.5) / 0.5
        const rotY = -(xCoord - 0.5) / 0.5
        const rotMagnitude = Math.max(Math.abs(rotX), Math.abs(rotY))
        const rotDeg = Math.round(rotationIntensity * rotMagnitude)
        transform = `perspective(1000px) rotate3d(${rotX}, ${rotY}, 0, ${rotDeg}deg)`

        shineX = xCoord
        shineY = yCoord
    }

    // - Shine
    const shineLeft = `calc(${100 - shineX * 100}% - ${shineSize / 2}px)`
    const shineTop = `calc(${100 - shineY * 100}% - ${shineSize / 2}px)`

    // - Thickness
    const layers = []
    for (let i = 0; i < depth; i++) {
        layers.push(
            <Box
                key={i}
                className={Classes.depthLayer}
                sx={{ transform: `translateZ(${i}px)` }}
            />
        );
    }

    return (
        <Box
            className={Classes.root}
            style={{
                transform,
            }}
            sx={{
                position: 'relative',
                transition: theme.transitions.create(
                    ['transform'],
                    {
                        duration: theme.transitions.duration.ease,
                    }
                ),
                '&, & > *': {
                    transformStyle: 'preserve-3d',
                },
                [`& > *, .${Classes.side}`]: {
                    borderRadius: '20px',
                    overflow: 'hidden',
                },
                [`.${Classes.depthLayer}`]: {
                    position: 'absolute',
                    top: -outlineSize,
                    left: -outlineSize,
                    right: -outlineSize,
                    bottom: -outlineSize,
                    backfaceVisibility: 'hidden',
                    bgcolor: '#FFFDEC',
                },
            }}
        >
            <Box
                sx={{
                    transform: `translateZ(${depth}px)`,
                    position: 'relative',
                }}
            >
                <Box
                    style={{
                        left: shineLeft,
                        top: shineTop,
                    }}
                    sx={{
                        position: 'absolute',
                        width: `${shineSize}px`,
                        height: `${shineSize}px`,
                        background: 'radial-gradient(circle, rgba(255,255,255,0.2) 0%, rgba(255,255,255,0.1) 25%, rgba(255,255,255,0) 50%)',
                        filter: 'blur(10px)',
                        borderRadius: '50%',
                        transition: theme.transitions.create(
                            ['top', 'left'],
                            {
                                duration: theme.transitions.duration.ease,
                            }
                        ),
                    }}
                />

                <Box ref={ref}>
                    <LibraryCardFront
                        user={user}
                    />
                </Box>
            </Box>
            {layers}
        </Box>
    )
}

function LibraryCardFront({
    user,
}: {
    user: User,
}) {
    const theme = useTheme()

    const [avatarLoaded, setAvatarLoaded] = React.useState(false)
    const avatarSideLengthPx = 24
    const avatarSideLength = `${24}px`

    return (
        <Box
            className={Classes.side}
            sx={{
                display: 'flex',
                flexDirection: 'column',
                width: '298px',
                color: 'white',
                borderRadius: '20px',
            }}
        >
            <Image
                src={user.imageUrl}
                alt={`Data Library image for ${user.name}`}
                style={{
                    width: '100%',
                    height: '250px',
                    objectFit: 'cover',
                }}
            />

            <Typography
                variant='h6SmallCaps'
                component='div'
                sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    gap: '10px',
                    px: '16px',
                    color: 'utilityWhite.main',
                    bgcolor: 'ground.grade10',
                }}
            >
                {user.role.toUpperCase()}

                <Box
                    sx={{
                        flexGrow: '1',
                        height: '1px',
                        bgcolor: 'utilityWhite.main',
                    }}
                />

                UCSF
            </Typography>

            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '16px',
                    p: '16px',
                    pt: '8px',
                    color: '#343434',
                    bgcolor: '#FFFDEC',
                }}
            >
                <Box>
                    <Typography variant='h5'>
                        {user.name}
                    </Typography>

                    <Typography variant='pMedium'>
                        {user.title}
                    </Typography>
                </Box>

                <Box
                    sx={{
                        height: '1px',
                        bgcolor: 'ground.grade10',
                        opacity: '0.2',
                    }}
                />

                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '16px',
                    }}
                >
                    {/* Avatar + email */}
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            gap: '6px',
                        }}
                    >
                        <Box
                            sx={{
                                minWidth: avatarSideLength,
                                maxWidth: avatarSideLength,
                                minHeight: avatarSideLength,
                                maxHeight: avatarSideLength,
                                borderRadius: '50%',
                                overflow: 'hidden',
                                position: 'relative',
                                '& > *': {
                                    width: '100%',
                                    minWidth: '100%',
                                    height: '100%',
                                    minHeight: '100%',
                                },
                            }}
                        >
                            <Box
                                sx={{
                                    display: 'flex',
                                    justifyContent: 'center',
                                    alignItems: 'center',
                                    bgcolor: user.color,
                                }}
                            >
                                <Typography
                                    variant='pBodyBoldWt'
                                    sx={{
                                        color: 'utilityHighlight.main',
                                    }}
                                >
                                    {user.name[0].toUpperCase()}
                                </Typography>
                            </Box>

                            {user.avatarUrl &&
                                <Image
                                    src={user.avatarUrl}
                                    alt={`Profile image for ${user.name}`}
                                    width={avatarSideLengthPx}
                                    height={avatarSideLengthPx}
                                    onLoad={() => setAvatarLoaded(true)}
                                    style={{
                                        opacity: avatarLoaded ? 1 : 0,
                                        position: 'absolute',
                                        top: 0,
                                        left: 0,
                                        transition: theme.transitions.create(
                                            ['opacity'],
                                            {
                                                easing: theme.transitions.easing.swell,
                                                duration: theme.transitions.duration.swell,
                                            },
                                        ),
                                    }}
                                />
                            }
                        </Box>

                        <Typography variant='pBody'>
                            {user.email}
                        </Typography>
                    </Box>

                    {/* Stats */}
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                        }}
                    >
                        {STATS.map(stat => (
                            <Box
                                key={stat.key}
                            >
                                <Typography
                                    variant='h6SmallCaps'
                                    component='div'
                                    sx={{
                                        px: '5px',
                                        color: 'ground.grade10',
                                    }}
                                >
                                    {stat.label.toUpperCase()}
                                </Typography>

                                <Typography
                                    variant='pMedium'
                                    component='div'
                                    sx={{
                                        px: '6px',
                                        py: '4px',
                                        color: '#343434'
                                    }}
                                >
                                    {user[stat.key] as number}
                                </Typography>
                            </Box>
                        ))}
                    </Box>

                    {/* Date + logo */}
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            justifyContent: 'space-between',
                            alignItems: 'center',
                        }}
                    >
                        <Typography
                            variant='pBodyMono'
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                width: 'fit-content',
                                gap: '9px',
                                px: '6px',
                                border: '1px solid #000',
                                borderRadius: '7px',
                                color: '#343434',
                            }}
                        >
                            <Box
                                sx={{
                                    py: '3px',
                                }}
                            >
                                DL
                            </Box>

                            <Box
                                sx={{
                                    width: '1px',
                                    bgcolor: 'ground.grade10',
                                }}
                            />

                            <Box
                                sx={{
                                    py: '3px',
                                }}
                            >
                                {user.joinDate.toLocaleDateString(undefined, {
                                    year: 'numeric',
                                    month: 'long',
                                    day: 'numeric',
                                })}
                            </Box>
                        </Typography>

                        <Image
                            src={logoDarkSrc}
                            alt='Data Library logo'
                            width={33}
                            height={33}
                        />
                    </Box>
                </Box>
            </Box>
        </Box>
    )
}

export default React.forwardRef(LibraryCard)