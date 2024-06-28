import { ThemeData } from "../themes/models";


export const ProjectStatus = {
    init: 'Initial Launch',
    assays: 'Assays',
    analysis: 'Analysis',
    public: 'Public Launch',
} as const

export const ProjectType = {
    coproject: 'CoProject',
    consortium: 'Consortium',
    copilot: 'CoPilot',
    collab: 'Collaboration',
    external: 'External',
} as const

export interface Project {
    name: string
    fullName: string
    heading?: string
    description: string
    fundingSource: string
    principalInvestigators: string[]
    status: typeof ProjectStatus[keyof typeof ProjectStatus]
    type: typeof ProjectType[keyof typeof ProjectType]
    species: string
    startDate: Date
    dataCollectionComplete: boolean
    userCount: number

    theme: ThemeData
}