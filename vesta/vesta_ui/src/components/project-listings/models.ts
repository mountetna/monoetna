import { ThemeData } from "../themes/models";


export enum ProjectStatus {
    team = 'Team',
    community = 'Community',
}

export enum ProjectType {
    coproject = 'CoProject',
    consortium = 'Consortium',
    copilot = 'CoPilot',
    collab = 'Collaboration',
    external = 'External',
    'hellman award' = 'Hellman Award',
}

export interface PrincipalInvestigator {
    name: string
    email: string
    title?: string
    imageUrl?: string
    profileUrl?: string
    color: string
    altColor: string
}

export interface Project {
    name: string
    fullName: string
    heading?: string
    description: string
    fundingSource?: string
    principalInvestigators: PrincipalInvestigator[]
    status: ProjectStatus
    type: ProjectType
    dataTypes: string[]
    sampleCount: number
    assayCount: number
    hasClinicalData: string,
    species: string
    startDate: Date
    dataCollectionComplete: boolean
    userCount: number
    theme: ThemeData
    href: string
}

export enum ExternalProjectStatus {
    init = 'Initial Launch',
    sampling = 'Sampling',
    analysis = 'Data Analysis',
    public = 'Community',
}

export function getExternalProjectStatus(proj: Project): ExternalProjectStatus {
    if (proj.status === ProjectStatus.community) {
        return ExternalProjectStatus.public
    } else if (proj.assayCount) {
        return ExternalProjectStatus.analysis
    } else if (proj.sampleCount) {
        return ExternalProjectStatus.sampling
    } else {
        return ExternalProjectStatus.init
    }
}

export enum ProjectHeadingInfoSet {
    default = 'Default',
    dataTypes = 'Data Types',
    pis = 'Principal Investigators',
}

export const PROJECTS_SEARCH_PARAMS_KEY = 'projects'

export interface ProjectsSearchParamsControls {
    viewSet?: string
    filterMethod?: string
    page?: number
}

export interface ProjectsSearchParamsState {
    filters?: Record<string, string[]>
    controls?: ProjectsSearchParamsControls
}