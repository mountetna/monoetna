'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material/styles';

import AboutItem, { Link } from './about-item'

interface AboutItemProps {
    header: string
    body: string
    link?: Link
    imageSrc: string
}

export default function AboutCarousel({ items }: { items: AboutItemProps[] }) {
    const theme = useTheme()

    const [carouselIndex, setCarouselIndex] = React.useState<number>(0)

    const handleChangeIndex = (indexNext: number, indexCurrent: number) => {
        setCarouselIndex(indexNext)
    }

    return (
        <Box
            sx={(theme) => ({

            })}
        >
            {/* <EnhancedSwipeableViews
                slideStyle={{ padding: '0 10px' }}
                containerStyle={{ padding: '0 10px' }}
                enableMouseEvents
                springConfig={{
                    duration: `${theme.transitions.duration.quint / 1000}s`,
                    easeFunction: theme.transitions.easing.quint,
                    delay: '0s',
                }}
                index={carouselIndex}
                onChangeIndex={handleChangeIndex}
            >
                {items.map((item) => {
                    return (
                        <React.Fragment key={item.header}>
                            <AboutItem
                                header={item.header}
                                body={item.body}
                                link={item.link}
                                imageSrc={item.imageSrc}
                            />
                        </React.Fragment>
                    )
                })}
            </EnhancedSwipeableViews> */}
        </Box>
    )
}