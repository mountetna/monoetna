'use client'

import * as React from 'react'
import Typography from '@mui/material/Typography';
import { styled, alpha, Box } from '@mui/system'
import TextInput from '../inputs/text-input';
import Toggle from '../inputs/toggle';
import Slider from '../inputs/slider';
import hideIcon from '/public/images/icons/hide.svg'
import showIcon from '/public/images/icons/show.svg'
import caretUpIcon from '/public/images/icons/caret-up.svg'
import caretDownIcon from '/public/images/icons/caret-down.svg'
import searchDarkIcon from '/public/images/icons/search.svg'
import infoSmallIcon from '/public/images/icons/info-small.svg'
import Autocomplete from '@/components/searchable-list/controls/autocomplete';
import ButtonBase from '@mui/material/ButtonBase';
import ProjectPI from './project-pi';
import FilterPill from '../searchable-list/filter-pill';
import FilterChip from '../searchable-list/filter-chip';
import Image from 'next/image';
import { FormControl } from '@mui/base';
import { ProjectExplorerContext } from '@/components/project-explorer/context';
import { DATA_TYPES } from '@/lib/fixtures';
import { projectDataTypes } from '@/lib/utils/filters';

const BasicToggle = ({on, labelOn, labelOff, altOn, altOff, iconOn, iconOff, onClick}:{
  on: boolean;
  iconOn: string;
  altOn: string;
  labelOn: string;
  iconOff: string;
  altOff: string;
  labelOff: string;
  onClick: Function;
}) => {
    return <ButtonBase
      onClick={onClick}
      aria-label={ on ? labelOn : labelOff }
      sx={{ aspectRatio: '1/1', padding: '6px' }}
    >
      <Image
          src={on ? iconOn : iconOff }
          alt={on ? altOn : altOff }
          width='26'
          height='26'
      />
    </ButtonBase>
}

const ShowHideToggle = ({shown, icon, label, onClick}:{
  shown: boolean;
  onClick: Function;
}) => {
    return <BasicToggle
      on={shown}
      onClick={onClick}

      labelOn={ 'Filter enabled' }
      iconOn={ showIcon }
      altOn={'Open eye'}

      labelOff={ 'Filter disabled' }
      iconOff={ hideIcon }
      altOff={'Closed eye'}
    />
}

const Pill = ({value,text}) => {
  return <Box sx={theme => ({
      width: '90px',
      display: 'flex',
      flexDirection: 'column',
      alignItems: 'center',
      '& input': {
        width: '100%',
        py: '5px',
        mb: '5px',
        textAlign: 'center',
        bgcolor: 'ground.grade100',
        ...theme.typography.pBodyMono
      }
    })} >
    <FormControl>
      <TextInput value={value} gradVariant='mid'/>
    </FormControl>
    <Typography variant='p3XSBoldWt' color='ground.grade50'>{text}</Typography>
  </Box>
}

const MinMaxSlider = ({min,max, value, setValue}) => {
  return <Box sx={{ display: 'flex' }} >
    <Pill value={value[0]} text="MIN"/>
    <Box sx={{ flex: '1 1 auto', px: '15px' }}>
      <Slider min={min} max={max} value={value} onChange={setValue} />
    </Box>
    <Pill value={value[1]} text="MAX"/>
  </Box>
}

const OpenCloseToggle = ({open, icon, label, onClick}:{
  open: boolean;
  onClick: Function;
}) => {
    return <BasicToggle
      on={open}
      onClick={onClick}

      labelOn={ 'Show filter options' }
      iconOn={ caretUpIcon }
      altOn={'Caret up'}

      labelOn={ 'Hide filter options' }
      iconOff={ caretDownIcon }
      altOff={'Caret down'}
    />
}

const Filter = ({title, highlight, children}:{
  title: string;
  children: React.ReactNode;
}) => {
  const [ active, setActive ] = React.useState(false);
  const [ fold, setFold ] = React.useState(true);

  const { state: { visibleColumns }, toggleColumnVisibility } = React.useContext(ProjectExplorerContext);

  const visible = visibleColumns.includes(title); 

  return <Box sx={{ py: '8px', my: '12px' }}>
    <Box sx={{ display: 'flex', alignItems: 'center'}} >
      { highlight && <Box sx={{
        bgcolor: 'green.grade75',
        width: '6px',
        height: '6px',
        marginRight: '5px',
        borderRadius: '50%',
      }}/> }
      <Box sx={{ flex: '1 1 auto' }}>
        <Typography variant="pBodyMediumWt">{title}</Typography>
      </Box>
      <Box sx={{ width: '76px' }}>
        <ShowHideToggle shown={visible} onClick={ () => toggleColumnVisibility(title) }/>
        <OpenCloseToggle open={!fold} onClick={ () => setFold(!fold) }/>
      </Box>
    </Box>
    {
      !fold && <Box sx={{ marginTop: '5px' }}>
        { children }
      </Box>
    }
  </Box>
}

const BasicFilter = ({title, collect, filter, render, id}) => {
  const { state: { projectData, filters, filterItemSet }, updateFilter, setFilterItems } = React.useContext(ProjectExplorerContext);

  const filterItems = filterItemSet[title] || [];
  const items = React.useMemo(
    () => {
      let items = {};
      projectData.forEach( project => collect(project,items) );
      return Object.values(items).sort( (a,b) => id(a).localeCompare(id(b)));
    }, [ projectData ]
  );

  const handleChangeFilterItems = (filterItems) => {
    updateFilter(title,
      !filterItems.length ? null :
      (project, matchAllFilters) => filterItems[
        matchAllFilters ? 'every' : 'some'
      ](
        filterItem => filter(filterItem, project, matchAllFilters)
      ),
      filterItems
    );
  };

  const handleClickRemoveFilterItem = React.useCallback(
    (filterItem) => {
      handleChangeFilterItems(filterItems.filter(f => f !== filterItem))
    }, [filterItems]
  );

  return <Filter title={title} highlight={ title in filters }>
    <Autocomplete
      size='small'
      multiple
      filterSelectedOptions
      icon={searchDarkIcon}
      options={items}
      onChange={(_, value, reason) => {
          if (reason === 'removeOption') return;
          handleChangeFilterItems(value);
      }}
      value={filterItems}
      getOptionLabel={ (option) => id(option) }
      getOptionKey={(option: FilterItem) => id(option)}
      isOptionEqualToValue={(option, value) => id(option) === id(value)}
      renderOption={ render }
      renderNoResults={() => (
          <Typography variant='pMedium'>
              No results
          </Typography>
      )}
      placeholder={title}/>
    {filterItems.length > 0 && <Box
        sx={{
            display: 'flex',
            flexDirection: 'row',
            flexWrap: 'wrap',
            columnGap: '10px',
            rowGap: '16px',
            pt: '16px'
        }}
    >
        {filterItems.map((item => (
            <FilterPill
                key={id(item)}
                label={id(item)}
                removeable
                onClickRemove={() => handleClickRemoveFilterItem(item)}
            />
        )))}
    </Box>}
  </Filter>
}

const InvestigatorsFilter = () => {
  return <BasicFilter
    title='Investigators'
    collect={
      (project, items) => {
        project.principalInvestigators.forEach(pi => items[pi.name] = pi);
      }
    }
    filter={
      (filterItem, project, matchAllFilters) => project.principalInvestigators.some(pi=> pi.name == filterItem.name)
    }
    render={
      params => (
        <ProjectPI
          data={params.option}
          showAvatar
          showNameAndTitle
          variant='filled'
        />
      )
    }
    id={ pi => pi.name }/>
}

const DEFAULT_RANGE=[0,500];
const DataTypesFilter = () => {
  return <BasicFilter
    title='Data types'
    collect={
        (project,items) => {
          const dataTypes = projectDataTypes(project);
          dataTypes.forEach( dt => items[dt] = dt );
        }
    }
    filter={
      (filterItem, project, matchAllFilters) => {
        return projectDataTypes(project).includes(filterItem)
      }
    }
    render={
      params => <Typography variant='pBodyMediumWt' key={params.option}>{params.option}</Typography>
    }
    id={ item => item }
  />
}

const ThemeFilter = () => {
  return <BasicFilter
    title='Theme'
    collect={
        (project,items) => {
          items[project.theme.name] = project.theme
        }
    }
    filter={
      (filterItem, project, matchAllFilters) => project.theme.name == filterItem.name
    }
    render={
      params => (
        <Typography variant='pMedium'>{params.option.name}</Typography>
      )
    }
    id={ item => item.name }
  />
}

const SUBJECTS_RANGE = [0,500];
const SubjectsFilter = () => {
  const { state: { projectData, filters, filterItemSet }, updateFilter } = React.useContext(ProjectExplorerContext);

  const subjectCounts = filterItemSet['Subjects'] || SUBJECTS_RANGE;

  const handleChangeSubjectCounts = (subjectCounts) => {
    updateFilter('Subjects',
      subjectCounts.every((e,i) => e == SUBJECTS_RANGE[i]) ? null :
      (project, matchAllFilters) => project.subjectCount >= subjectCounts[0] && (project.subjectCount <= subjectCounts[1] || subjectCounts[1] == SUBJECTS_RANGE[1]),
      subjectCounts
    );
  }

  return <Filter title="Subjects" highlight={ 'Subjects' in filters } >
    <MinMaxSlider min={SUBJECTS_RANGE[0]} max={SUBJECTS_RANGE[1]} value={subjectCounts} setValue={
      (e,value) => handleChangeSubjectCounts(value)
    }/>
  </Filter>
}

const TypeFilter = () => {
  return <BasicFilter
    title='Type'
    collect={
        (project,items) => {
          items[project.type] = project.type
        }
    }
    filter={
      (filterItem, project, matchAllFilters) => project.type == filterItem
    }
    render={
      params => (
        <Typography variant='pMedium'>{params.option}</Typography>
      )
    }
    id={ item => item }
  />
}

const PhaseFilter = () => {
  return <BasicFilter
    title='Phase'
    collect={
        (project,items) => {
          items[project.status] = project.status
        }
    }
    filter={
      (filterItem, project, matchAllFilters) => project.status == filterItem
    }
    render={
      params => (
        <Typography variant='pMedium'>{params.option}</Typography>
      )
    }
    id={ item => item }
  />
}

const currentYear = new Date().getFullYear();
const YEAR_RANGE = [2020,currentYear];

const YearFilter = () => {
  const { state: { projectData, filters, filterItemSet }, updateFilter } = React.useContext(ProjectExplorerContext);

  const dateRange = filterItemSet['Year'] || YEAR_RANGE;

  const handleChangeDateRange = (dateRange) => {
    updateFilter('Year',
      dateRange.every((e,i) => e == YEAR_RANGE[i]) ? null :
      (project, matchAllFilters) => project.startDate.getFullYear() >= dateRange[0] && project.startDate.getFullYear() <= dateRange[1],
      dateRange
    );
  }

  return <Filter title="Year" highlight={ 'Year' in filters } >
    <MinMaxSlider min={YEAR_RANGE[0]} max={YEAR_RANGE[1]} value={dateRange} setValue={
      (e,value) => handleChangeDateRange(value)
    }/>
  </Filter>
}

const FilterList = () => {
  const { state: { matchAllFilters }, setMatchAllFilters } = React.useContext(ProjectExplorerContext);
  
  return <>
    <Box sx={{
      display: 'flex',
      alignItems: 'center',
      pr: '5px',
      gap: '5px'
    }}>
      <Typography variant='pBodyMediumWt'>Match all filters</Typography>
      <Image
          src={infoSmallIcon}
          alt={'Info'}
          width='18'
          height='18'
      />
      <Box sx={{ flex: '1 1 auto' }}/>
      <Toggle active={ matchAllFilters } setActive={ setMatchAllFilters }/>
    </Box>
    <InvestigatorsFilter/>
    <DataTypesFilter/>
    <ThemeFilter/>
    <SubjectsFilter/>
    <TypeFilter/>
    <PhaseFilter/>
    <YearFilter/>
  </>
}

export default FilterList;
