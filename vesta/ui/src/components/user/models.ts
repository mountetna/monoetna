import { StaticImageData } from "next/image"


export interface User {
    name: string
    email: string
    title: string
    role: string
    imageUrl: StaticImageData
    avatarUrl?: string
    color: string
    joinDate: Date,
    contributions: number
    projectMemberships: number
}