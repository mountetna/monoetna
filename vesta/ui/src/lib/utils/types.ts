import * as React from 'react'
import { TypographyPropsVariantOverrides } from '@mui/material/Typography';
import { OverridableStringUnion } from '@mui/types'
import { Variant } from '@mui/material/styles/createTypography';


export type ValueOf<T> = Required<T>[keyof T];

export function forwardRef<T, P = {}>(
  render: (props: P, ref: React.Ref<T>) => React.ReactNode
): (props: P & React.RefAttributes<T>) => React.ReactNode {
  // @ts-ignore
  return React.forwardRef(render) as any;
}

export type TypographyVariant = OverridableStringUnion<"inherit" | Variant, TypographyPropsVariantOverrides>