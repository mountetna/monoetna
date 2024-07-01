import { StaticImageData } from "next/image"


interface ThemeImageComponents {
    filtered: StaticImageData
    projectBackground: StaticImageData
}


export interface ThemeData {
    name: string
    description: string
    projectCount: number
    projectsLink: string
    color: string
    textColor: 'light' | 'dark'
    imageComponents: ThemeImageComponents
    coverImage: StaticImageData
    icon: StaticImageData
}