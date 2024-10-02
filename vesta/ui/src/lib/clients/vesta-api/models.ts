import { StatsTimeseries } from "@/components/stats/library-stats"

export interface ApiStatsInstance {
    id: string
    user_count: number
    file_count: number
    byte_count: number
    subject_count: number
    sample_count: number
    assay_count: number
    recorded_at: Date
}

export interface ApiProjectStats extends ApiStatsInstance {
    name: string
    clinical_data_count: number
}

export interface Stats {
    global: StatsTimeseries
    byProjectName: Record<string, StatsTimeseries>
}

interface PrincipalInvestigator {
    name: string
    email: string
    title?: string
    image_url?: string
    profile_url?: string
}

export interface ApiProjectInfo {
    id: string
    name: string
    full_name: string
    description: string
    funding_source?: string
    principal_investigators: PrincipalInvestigator[]
    status: string
    type: string
    data_types: string[]
    species: string
    start_date: Date  // really just the year is valid
    theme: string  // just the name
    data_collection_complete: boolean
    created_at: Date
    updated_at: Date
}