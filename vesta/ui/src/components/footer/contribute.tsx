'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import Image from 'next/image';
import { useTheme } from '@mui/material';

import contributeImage from '/public/images/footer/contribute.png'
import TextInput from '../inputs/text-input';
import ArrowLinkButton from '../inputs/arrow-link-button';
import { FormControl } from '@mui/base';


export default function Contribute({ }: {}) {
    const theme = useTheme()

    const [inputVal, setInputVal] = React.useState<string>('')

    const formRef = React.useRef<HTMLFormElement>()
    const handleSubmitForm = (event?: React.FormEvent) => {
        const formEl = formRef.current

        if (!formEl) return
        if (!formEl.reportValidity()) {
            if (!event) {
                return
            }
        }

        event && event.preventDefault()

        console.log('here', inputVal)
    }

    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'column',
                gap: '12px',
                [theme.breakpoints.up('tablet')]: {
                    flexDirection: 'row',
                    gap: '16px',
                },
                '& > *': {
                    height: '237.75px',
                    borderRadius: '30px',
                    overflow: 'hidden',
                    [theme.breakpoints.up('tablet')]: {
                        width: '50%',
                        height: '487.5px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        height: '485.5px',
                    },
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
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '30px',
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
                    ref={formRef}
                    component='form'
                    onSubmit={handleSubmitForm}
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        gap: '10px',
                        '& > *:first-child': {
                            display: 'flex',
                            flexGrow: '1',
                            [theme.breakpoints.up('desktop')]: {
                                minWidth: '393px',
                                flexGrow: 'unset',
                            },
                            '& > *': {
                                flexGrow: '1',
                            },
                            '& input': {
                                width: '100%',
                            },
                        },
                        '& > *:last-child': {
                            aspectRatio: '1',
                            height: '100%',
                            '& img': {
                                width: '30px',
                                height: '30px',
                            },
                        },
                    }}
                >
                    <FormControl>
                        <TextInput
                            value={inputVal}
                            onChange={(event) => setInputVal(event.currentTarget.value)}
                            placeholder='name@email.com'
                            gradeVariant='light'
                            required={true}
                            type='email'
                        />
                    </FormControl>

                    <ArrowLinkButton
                        onClick={() => handleSubmitForm()}
                    />
                </Box>
            </Box>
        </Box>
    )
}