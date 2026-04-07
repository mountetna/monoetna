'use client'

import * as React from 'react';
import Image from 'next/image';
import InputBase from '@mui/material/InputBase';
import Box from '@mui/system/Box';
import Typography from '@mui/material/Typography';
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';
import publicationsHeroImg from '/public/images/publications-hero-img.svg'
import searchDarkIcon from '/public/images/icons/search.svg'
import Pagination, { PaginationClasses } from '@/components/searchable-list/controls/pagination'
import { SxProps, useTheme } from '@mui/material';
import arrowUpRightDark from '/public/images/icons/arrow-up-right-dark.svg';
import Link from './link/link';
import { Publication } from '@/lib/clients/vesta-api/models';

export const PublicationsHero = () => {
  return (
    <Box sx={{
      mx: '100px',
      pt: '100px',
      height: '540px',
      background: `url(${publicationsHeroImg.src})`,
      backgroundPosition: '-100px -100px' }}>
      <Typography variant='h2'>Publication resources</Typography>
      <Box sx={{ mt: '32px', width: 650 }}>
        <Typography variant='pMedium'>
        The Data Library hosts a series of research projects containing multi-omic data and de-identified clinical data. Project data is shared in our research community to best leverage these rich datasets toward biological discovery.
        </Typography>
      </Box>
    </Box>
  );
}

const PubInfo = ({label, info}:{
  label: string;
  info: string|number;
}) => {
  const theme = useTheme();
  return <Box sx={{
    gap: '10px',
    display: 'flex'
  }}>
    <Typography variant='pBodyMediumWt' sx={{ color: theme.palette.ground.grade50 }}>{label}</Typography>
    <Typography variant='pBody' sx={{ color: theme.palette.ground.grade10 }}>{info}</Typography>
  </Box>
}

export const PublicationsTable = ({publications}:{
  publications: Publication[]
}) => {
  const theme = useTheme();
  const [ filterText, setFilterText ] = React.useState('');

  const filterTerms = filterText.split(' ');
  const filteredPublications = publications.filter(
    (pub:Publication) => {
      const pubText = `${pub.title} ${pub.authors} ${pub.publication_year} ${pub.project} ${pub.journal}`;
      return filterTerms.every(
        filterTerm => pubText.search(new RegExp(filterTerm, "i")) != -1
      )
    }
  ).sort( (pub1,pub2) => pub2.publication_year - pub1.publication_year );
  return (
    <Box sx={{ mx: '100px', mb: '100px' }}>
      <Typography variant='h2'>Publications</Typography>
      <Box sx={{ mt: '32px' }}>
        <Box sx={{ display: 'flex', justifyContent: 'space-between', borderBottom: `1px solid ${theme.palette.ground.grade50}`, pb: '17px', mb: '16px' }}>
          <Box
            sx={{
              width: '400px',
              height: '100%',
              overflowX: 'scroll',
              overflowY: 'clip',
              display: 'flex',
              gap: '10px',
              alignItems: 'right',
              bgcolor: 'utilityWhite.main',
              borderRadius: '30px',
              p: '7px 12px',
              border: '1px solid transparent',
              '&:focus-visible, &.Mui-focused': {
                border: `1px solid ${theme.palette.ground.grade75}`,
              },
              '&:hover, &:focus-visible, &.Mui-focused': {
                bgcolor: 'ground.grade100',
              },
              transition: theme.transitions.create(
                'all',
                {
                  easing: theme.transitions.easing.ease,
                  duration: theme.transitions.duration.ease,
                }
              ),
            }}
          >
            <Image
              src={searchDarkIcon}
              alt='Search icon'
            />

            <InputBase
              placeholder={'Search e.g. "Nature"'}
              value={filterText}
              onChange={ e => setFilterText(e.target.value) }
              sx={{
                border: 'none',
                outline: 'none',
                p: '0',
                background: 'inherit',
                ...theme.typography.pBody,
                color: 'ground.grade10',
                '&::placeholder': {
                  color: '#777777'
                },
              }}
            />
          </Box>
          <Pagination
              currentPage={0}
              onClickPrev={ () => {} }
              onClickNext={ () => {} }
              pageSize={filteredPublications.length}
              listSize={publications.length}
              listItemLabel='publications'
              showArrows={false}
          />
        </Box>
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: '16px'}}>
              {filteredPublications.map((pub:Publication,i) => (
                <Box
                  key={i}
                  sx={ { 
                    background: theme.palette.utilityWhite.main,
                    borderRadius: '16px',
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '24px',
                    p: '36px'
                  }}
                >
                  <Typography variant='h5BoldWt' sx={{ color: theme.palette.teal.grade75}}>{pub.title}</Typography>

                  <Box sx={{ display: 'flex', gap: '27px' }}>
                    <PubInfo label='Year' info={pub.publication_year}/>
                    <PubInfo label='Project' info={pub.project}/>
                    <PubInfo label='Journal' info={pub.journal}/>
                    <Link style={{ display: 'flex' }} target='_blank' href={pub.doi}>
                      <Typography variant='pBody' sx={{ color: theme.palette.ground.grade10 }}>Open link</Typography>
                      <Image
                        style={{ margin: '0px'}}
                        width={20}
                        height={20}
                        src={arrowUpRightDark}
                        alt='Arrow pointing up-right'
                      />
                    </Link>
                  </Box>

                  <Box sx={{ display: 'flex' }}>
                    <Box sx={{ flex: '0 0 100px' }}><Typography variant='pBodyMediumWt' sx={{ color: theme.palette.ground.grade50 }}>Authors</Typography></Box>
                    <Box sx={{ flex: '1 1 auto' }}>
                      <Typography variant='pBody' sx={{ color: theme.palette.ground.grade10 }}>{pub.authors}</Typography>
                    </Box>
                  </Box>
                </Box>
              ))}
        </Box>
      </Box>
    </Box>
  );
}
