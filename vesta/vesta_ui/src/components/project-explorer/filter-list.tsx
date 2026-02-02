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

const BasicFilter = ({title, filter, items, render, id}) => {
  const {
    state: { projectData, filters, filterItemSet },
    createFilter,
    updateFilterItems
  } = React.useContext(ProjectExplorerContext);

  // setup the filter in advance
  React.useEffect(
    () => {
      createFilter(title, filter, items);
    }, []
  )

  const filterItems = filterItemSet[title] || [];
  const projectItems = React.useMemo(
    () => {
      let projectItems = {};
      projectData.forEach( project => projectItems = { ...projectItems, ...items(project) } );
      return Object.values(projectItems).sort( (a,b) => id(a).localeCompare(id(b)));
    }, [ projectData ]
  );

  const handleChangeFilterItems = (filterItems) => {
    updateFilterItems(title, filterItems);
  };

  const handleClickRemoveFilterItem = React.useCallback(
    (filterItem) => {
      handleChangeFilterItems(filterItems.filter(f => f !== filterItem))
    }, [filterItemSet]
  );

  return <Filter title={title} highlight={ title in filterItemSet }>
    <Autocomplete
      size='small'
      multiple
      filterSelectedOptions
      icon={searchDarkIcon}
      options={projectItems}
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

const FreeFilter = () => {
  const {
    state: { projectData, filters, filterItemSet },
    createFilter,
    updateFilterItems
  } = React.useContext(ProjectExplorerContext);

  const filterItems = filterItemSet['free'];

  const handleClickRemoveFilterItem = React.useCallback(
    (filterItem) => {
      updateFilterItems('free', filterItems.filter(f => f !== filterItem));
    }, [filterItemSet]
  );

  if (!('free' in filterItemSet)) return null;

  return <Box>
    <Box sx={{ py: '10px' }}><Typography variant="pBodyMediumWt">Search keywords</Typography></Box>
    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: '10px'}}>
    {
      filterItems.map((item => (
            <FilterPill
                key={item}
                label={item}
                removeable
                onClickRemove={() => handleClickRemoveFilterItem(item)}
            />
        )))
    }
    </Box>
  </Box>
}

const InvestigatorsFilter = () => {
  return <BasicFilter
    title='Investigators'
    items={
      project => Object.fromEntries( project.principalInvestigators.map( pi => [ pi.name, pi ] ) )
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
    items={
      project => Object.fromEntries( projectDataTypes(project).map( dt => [ dt, dt ] ) )
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
    items={
      project => ({ [project.theme.name]: project.theme })
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
  return <RangeFilter
    title='Subjects'
    defaultRange={SUBJECTS_RANGE}
    rangeValue={p=>p.subjectCount}/>
}

const RangeFilter = ({title, defaultRange, rangeValue}) => {
  const { state: { projectData, filters, filterItemSet }, createFilter, updateFilterItems } = React.useContext(ProjectExplorerContext);

  const filter = (range, project, matchAllFilters) => {
    if (range.every((e,i) => e == defaultRange[i])) return true;

    return (rangeValue(project) >= range[0]
      && (rangeValue(project) <= range[1]
      || range[1] == defaultRange[1]))
  }

  React.useEffect( () => {
    createFilter(title, filter, null)
  }, [] );

  const handleChangeRange = (range) => {
    if (range.every((e,i) => e == defaultRange[i])) 
    updateFilterItems(title, null);
    else
    updateFilterItems(title, [ range ]);
  }

  const range = title in filterItemSet ? filterItemSet[title][0] : defaultRange;

  return <Filter title={title} highlight={ title in filterItemSet } >
    <MinMaxSlider min={defaultRange[0]} max={defaultRange[1]} value={range} setValue={
      (e,value) => handleChangeRange(value)
    }/>
  </Filter>
}

const TypeFilter = () => {
  return <BasicFilter
    title='Type'
    items={
      project => ({ [project.type]: project.type })
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
    items={
      project => ({ [project.status]: project.status })
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
  return <RangeFilter
    title='Year'
    defaultRange={YEAR_RANGE}
    rangeValue={p=>p.startDate.getFullYear()}/>
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
    <FreeFilter/>
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
