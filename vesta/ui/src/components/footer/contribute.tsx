'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import Image from 'next/image';

import contributeImage from '/public/images/footer/contribute.png'
import TextInput from '../inputs/text-input';


export default function Contribute({ }: {}) {
    const [inputVal, setInputVal] = React.useState<string>('')

    const handleSubmitForm = (event: React.FormEvent) => {
        event.preventDefault()

        console.log('here', inputVal)
    }

    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'column',
                gap: '12px',
                '& > *': {
                    // width: '100%',
                    // height: '237.75px',
                    aspectRatio: '377 / 237.75',
                    borderRadius: '30px',
                    overflow: 'hidden',
                },
            }}
        >
            <Box>
                <Image
                    src={contributeImage}
                    alt='Spatial transcriptomic image analysis for Tonsil project'
                    style={{
                        width: '100%',
                        height: '100%',
                        objectFit: 'cover',
                    }}
                />
            </Box>

            <Box
                sx={{
                    p: '16px',
                    bgcolor: '#3A7FA3',
                }}
            >
                <Typography
                    variant='h4'
                    component='div'
                    sx={{
                        color: 'utilityWhite.main'
                    }}
                >
                    Researching something groundbreaking? Add your data to the library.
                </Typography>

                <Box
                    component='form'
                    onSubmit={handleSubmitForm}
                >
                    <TextInput
                        value={inputVal}
                        onChange={(event) => setInputVal(event.currentTarget.value)}
                        placeholder='name@email.com'
                    />
                </Box>
            </Box>
        </Box>
    )
}