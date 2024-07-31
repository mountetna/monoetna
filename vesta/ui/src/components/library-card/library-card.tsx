import * as React from 'react';
import Box from '@mui/system/Box'
import { Typography, useTheme } from '@mui/material';
import Image from 'next/image';

import { User } from '../user/models';

import logoDarkSrc from '/public/images/logo/logo-dark.svg'


export enum Classes {
    root = 'library-card'
}


interface Stat {
    key: keyof User
    label: string
}

const STATS: Stat[] = [
    { key: 'contributions', label: 'Contributions', },
    { key: 'projectMemberships', label: 'Projects', },
]

interface Props {
    user: User,
}

function LibraryCard(props: Props, ref?: React.ForwardedRef<HTMLElement>) {
    const theme = useTheme()
    const user = props.user

    const [avatarLoaded, setAvatarLoaded] = React.useState(false)
    const avatarSideLengthPx = 24
    const avatarSideLength = `${24}px`

    return (
        <Box
            ref={ref}
            className={Classes.root}
            sx={{
                display: 'flex',
                flexDirection: 'column',
                width: '298px',
                borderRadius: '20px',
                overflow: 'hidden',
                color: 'white',
            }}
        >
            <Box
                sx={{
                    border: '1px solid white',
                    height: '200px',
                }}
            >
                image placeholder
            </Box>

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
                {user.role}

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

                            {user.imageUrl &&
                                <Image
                                    src={user.imageUrl}
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
                            <Box>
                                <Typography
                                    variant='h6SmallCaps'
                                    component='div'
                                    sx={{
                                        px: '5px',
                                        color: 'ground.grade10',
                                    }}
                                >
                                    {stat.label}
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
                                {user.joinDate.toLocaleDateString()}
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