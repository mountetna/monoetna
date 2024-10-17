import { defaultDict } from "@/lib/utils/object"
import { ApiProjectInfo, ApiProjectStats, Stats, ApiStatsInstance, SendContactStatus } from "./models"
import { StatsTimeseries } from "@/components/stats/library-stats"

export class VestaApiClient {
    baseUrl: string

    constructor() {
        this.baseUrl = process.env.API_URL
    }

    async fetchStats(): Promise<Stats> {
        const [
            globalStatsRes,
            projectStatsRes,
        ] = await Promise.all([
            fetch(this.#createUrl('stats')),
            fetch(this.#createUrl('stats/projects')),
        ])

        let globalStats: ApiStatsInstance[]
        let projectStats: ApiProjectStats[]
        [
            globalStats,
            projectStats,
        ] = await Promise.all([
            globalStatsRes.json(),
            projectStatsRes.json(),
        ])

        const globalTimeseries = this.#createEmptyTimeseries()

        globalStats.sort((a, b) => {
            const [aDate, bDate] = [new Date(a.recorded_at), new Date(b.recorded_at)]
            return aDate.getTime() - bDate.getTime()
        })

        for (const stats of globalStats) {
            const date = new Date(stats.recorded_at)

            for (const [k, v] of (Object.entries(stats) as [keyof ApiStatsInstance, number][])) {
                if (k == 'id' || k == 'recorded_at') {
                    continue
                }
                globalTimeseries[this.#apiStatsKeyToTimeseriesKey(k)].push({
                    value: v,
                    date,
                })
            }
        }

        projectStats.sort((a, b) => {
            const [aDate, bDate] = [new Date(a.recorded_at), new Date(b.recorded_at)]
            return aDate.getTime() - bDate.getTime()
        })

        const byProjectName = defaultDict<string, StatsTimeseries>(_ => this.#createEmptyTimeseries())

        for (const stats of projectStats) {
            const date = new Date(stats.recorded_at)
            const name = stats.name

            for (const [k, v] of (Object.entries(stats) as [keyof ApiProjectStats, number][])) {
                if (k == 'id' || k == 'recorded_at' || k == 'name' || k == 'clinical_data_count') {
                    continue
                }
                byProjectName[name][this.#apiStatsKeyToTimeseriesKey(k)].push({
                    value: v,
                    date,
                })
            }
        }



        return {
            global: globalTimeseries,
            byProjectName,
        }
    }

    async fetchProjects(): Promise<ApiProjectInfo[]> {
        const res = await fetch(this.#createUrl('projects'))
        const apiProjects: ApiProjectInfo[] = await res.json()

        return apiProjects.map(proj => ({
            ...proj,
            start_date: new Date(proj.start_date)
        }))
    }

    #createUrl(path: string): string {
        return (new URL(path, this.baseUrl)).href
    }

    #createEmptyTimeseries(): StatsTimeseries {
        return {
            bytes: [],
            assays: [],
            subjects: [],
            files: [],
            samples: [],
            users: [],
        }
    }

    #apiStatsKeyToTimeseriesKey(apiKey: keyof ApiStatsInstance): keyof StatsTimeseries {
        // @ts-ignore
        return apiKey.split('_')[0] + 's'
    }

    async sendContributeEmail(email: string): Promise<SendContactStatus> {
        let res: Response
        try {
            res = await fetch(this.#createUrl('contact'), {
                method: "POST",
                body: JSON.stringify({ requester_email: email }),
            })
        } catch (e) {
            return {
                status: 'error',
                message: String(e),
            }
        }

        if (res && res.status === 200) {
            return {
                status: 'success',
            }
        }

        return {
            status: 'error',
            message: 'Unknown error'
        }
    }
}