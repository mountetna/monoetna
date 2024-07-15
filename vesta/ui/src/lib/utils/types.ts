import * as React from 'react'


export type ValueOf<T> = Required<T>[keyof T];

export function forwardRef<T, P = {}>(
    render: (props: P, ref: React.Ref<T>) => React.ReactNode
  ): (props: P & React.RefAttributes<T>) => React.ReactNode {
    return React.forwardRef(render) as any;
  }