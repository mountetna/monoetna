import { VestaApiClient } from '@/lib/clients/vesta-api/client';
import { defaultDict } from '@/lib/utils/object';
import { Project, ProjectStatus, ProjectType } from '@/components/project-explorer/models';
import { VIDEOS, ABOUT_ITEMS, THEMES } from '@/lib/fixtures';
import { faker } from '@faker-js/faker'
import _ from 'lodash';
import { Stats } from '@/components/stats/stats-carousel'

export async function getData() {
  const apiClient = new VestaApiClient()

  const [
    stats,
    apiProjects,
  ] = await Promise.all([
    apiClient.fetchStats(),
    apiClient.fetchProjects(),
  ])

  const themeProjectsCount = defaultDict<string, number>(_ => 0)

  const projects: (Project | undefined)[] = apiProjects.map(proj => {
    const status = Object.values(ProjectStatus).find(k => k.toUpperCase() === proj.status.toUpperCase())
    const type = Object.values(ProjectType).find(k => k.toUpperCase() === proj.type.toUpperCase())
    const theme = THEMES.find(theme => theme.name.toUpperCase() === proj.theme.toUpperCase())
    if (!status || !type || !theme) {
      return
    }

    themeProjectsCount[theme.name.toUpperCase()] += 1

    const dataTypes = proj.data_types.filter(dt => dt.toUpperCase() !== 'project'.toUpperCase())

    const latestCount = (stat:keyof Stats):number => {
      if (!(proj.name in stats.byProjectName)) return 0;

      const projStats = stats.byProjectName[proj.name];

      if (!(stat in projStats)) return 0;

      return projStats[stat].at(-1)?.value || 0;
    }

    return {
      name: proj.name,
      fullName: _.words(proj.full_name).join(' '),
      description: proj.description,
      fundingSource: proj.funding_source,
      principalInvestigators: proj.principal_investigators.map(pi => {
        const theme = faker.helpers.arrayElement(THEMES)

        return {
          name: pi.name,
          email: pi.email,
          title: pi.title,
          imageUrl: pi.photo_url,
          profileUrl: pi.profile_url,
          color: theme.color,
          altColor: theme.textColor === 'light' ? 'utilityHighlight.main' : 'ground.grade10',
        }
      }),
      status,
      type,
      dataTypes,
      sampleCount: latestCount('samples'),
      subjectCount: latestCount('subjects'),
      assayCount: latestCount('assays'),
      hasClinicalData: latestCount('assays') > 0 ? 'Yes' : 'No',
      species: proj.species,
      startDate: proj.start_date,
      dataCollectionComplete: proj.data_collection_complete,
      userCount: latestCount('users'),
      theme,
      href: (new URL(`/${proj.name}`, process.env.TIMUR_URL)).href
    }
  })

  const themes = THEMES.map(theme => ({
    ...theme,
    projectCount: themeProjectsCount[theme.name.toUpperCase()],
  }))

  return {
    heroVideos: VIDEOS,
    stats: stats.global,
    aboutItems: ABOUT_ITEMS,
    themes,
    projects: (projects.filter(p => p) as Project[]),
    accessUrl: process.env.TIMUR_URL
  }
}
