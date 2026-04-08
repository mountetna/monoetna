import * as React from 'react';
import { PublicationsTable, PublicationsHero } from '@/components/publications';
import { getData } from '@/lib/clients/vesta-api/request';

export default async function Publications() {
  const data = await getData();

  const publications = data.projects.map(project => project.publications).flat().filter(_=>_);

  return <>
    <PublicationsHero/>
    <PublicationsTable publications={publications}/>
  </>
}
