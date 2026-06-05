import * as React from 'react'
import { Box, Fade, useTheme } from '@mui/material';
import { Group } from '@visx/group';
import { Pack, hierarchy } from '@visx/hierarchy';
import { useParentSize } from '@visx/responsive'
import _ from 'lodash';


interface UserCountItem {
    radius: number
}


export default function ProjectUserCountPack({
    userCount,
    color,
    visible,
}: {
    userCount: number,
    color: string,
    visible: boolean,
}) {
    const theme = useTheme()

    let {
        parentRef: userCountPackContainerRef,
        width: userCountPackContainerWidth,
        height: userCountPackContainerHeight,
    } = useParentSize({ debounceTime: 500, })

    // Add some variability to similar-looking packs
    const packMarginPx = React.useMemo(() => _.random(3, 15, false), [])
    const rotationDeg = React.useMemo(() => _.random(0, 360, false), [])

    userCountPackContainerWidth -= packMarginPx * 2
    userCountPackContainerHeight -= packMarginPx * 2

    // Make pack item size dependent on container width
    // so we can consistently size the center gap as a portion of container width
    const packItemAreaPx = Math.PI * ((userCountPackContainerWidth / 2) ** 2) / userCount
    const packItemRadiusPx = Math.sqrt(packItemAreaPx) / Math.PI
    const gapItemRadius = 0.33 * userCountPackContainerWidth / 2

    const packItems: UserCountItem[] = Array(userCount).fill({ radius: packItemRadiusPx })
    packItems.unshift({ radius: gapItemRadius })

    const pack = {
        children: packItems,
        name: 'root',
        radius: 0,
        distance: 0,
    }

    const packRoot = hierarchy<UserCountItem>(pack)
        .sum(d => d.radius ** 2)

    return (
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
            {userCountPackContainerWidth > 10 && (
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
                                <Group
                                    style={{
                                        position: 'relative',
                                    }}
                                >
                                    {circles.map((circle, i) => (
                                        <Fade
                                            key={`fade-${i}`}
                                            in={visible}
                                            easing={theme.transitions.easing.ease}
                                            timeout={theme.transitions.duration.ease}
                                            style={{
                                                transitionDelay: `${theme.transitions.duration.ease / circles.length * i + theme.transitions.duration.ease / 1.1}ms`,
                                            }}
                                        >
                                            <circle
                                                key={`circle-${i}`}
                                                r={circle.r}
                                                cx={circle.x}
                                                cy={circle.y}
                                                fill={color}
                                            />
                                        </Fade>
                                    ))}
                                </Group>
                            )
                        }}
                    </Pack>
                </svg>
            )}
        </Box>
    )
}