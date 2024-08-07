'use client'

import * as React from 'react';
import { Input as BaseInput, InputProps } from '@mui/base/Input';
import Typography from '@mui/material/Typography';
import { styled } from '@mui/material';
import _ from 'lodash'

// import { styled } from '@/lib/utils/types';


type GradeVariant = 'light' | 'mid'

const BgColorDefault: Record<GradeVariant, string> = {
    light: 'utilityWhite.main',
    mid: 'utilityHighlight.main',
}

const BgColorHover: Record<GradeVariant, string> = {
    light: 'ground.grade100',
    mid: 'ground.grade100',
}

const BgColorFocus: Record<GradeVariant, string> = {
    light: 'utilityWhite.main',
    mid: 'utilityHighlight.main',
}

interface _Props {
    gradeVariant: GradeVariant
}

type Props = InputProps & _Props

const Input = React.forwardRef(function CustomInput(
    props: Props,
    ref: React.ForwardedRef<HTMLDivElement>,
) {
    const gradeVariant = props.gradeVariant
    const inputProps = _.omit(props, ['gradeVariant'])

    const InputElement = React.useMemo(() => styled('input')(
        ({ theme }) => theme.unstable_sx({
            px: '16px',
            py: '8px',
            border: '1px solid transparent',
            borderRadius: '30px',
            color: 'ground.grade10',
            bgcolor: BgColorDefault[gradeVariant],
            ...theme.typography.pLarge,
            transition: theme.transitions.create(
                'all',
                {
                    easing: theme.transitions.easing.ease,
                    duration: theme.transitions.duration.ease,
                }
            ),

            [theme.breakpoints.up('tablet')]: {
                px: '16px',
                py: '14px',
            },

            '&::placeholder': {
                color: '#777777',
            },

            '&:hover': {
                bgcolor: BgColorHover[gradeVariant],
            },

            '&:focus': {
                borderColor: 'ground.grade50',
                bgcolor: BgColorFocus[gradeVariant],
            },

            // firefox
            '&:focus-visible': {
                outline: 0,
            },
        })
    ), [gradeVariant])

    return (
        // @ts-ignore
        <BaseInput slots={{ input: InputElement }} {...inputProps} ref={ref} />
    )
});

export default function TextInput(props: Props) {
    return (
        <Input {...props} />
    )
}