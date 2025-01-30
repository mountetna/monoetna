'use client'

import * as React from 'react'

import { User } from './models'


export const UserContext = React.createContext<User | null>(null)

export const useUser = () => React.useContext(UserContext)

export function UserContextProvider({
    value,
    children,
}: {
    value: User | null,
    children: React.ReactNode,
}) {
    return (
        <UserContext.Provider value={value}>
            {children}
        </UserContext.Provider>
    )
}