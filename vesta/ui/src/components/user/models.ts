export interface User {
    name: string
    email: string
    title: string
    role: string
    imageUrl?: string
    color: string
    joinDate: Date,
    contributions: number
    projectMemberships: number
}