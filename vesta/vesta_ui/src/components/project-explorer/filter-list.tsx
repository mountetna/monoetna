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
import Image from 'next/image';
import { FormControl } from '@mui/base';

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
      ariaLabel={ on ? labelOn : labelOff }
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

const MinMaxSlider = ({min,max, currMin, currMax}) => {
  const minPercent = (currMin - min) / (max - min);
  const maxPercent = (currMax - min) / (max - min);
  return <Box sx={{ display: 'flex' }} >
    <Pill value={currMin} text="MIN"/>
    <Box sx={{ flex: '1 1 auto', px: '15px' }}>
      <Slider min={min} max={max} value={[currMin, currMax]} />
    </Box>
    <Pill value={currMax} text="MAX"/>
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

const Filter = ({title, children}:{
  title: string;
  children: React.ReactNode;
}) => {
  const [ active, setActive ] = React.useState(false);
  const [ shown, setShown ] = React.useState(false);

  return <Box sx={{ py: '8px', my: '12px' }}>
    <Box sx={{ display: 'flex', alignItems: 'center'}} >
      <Box sx={{ flex: '1 1 auto' }}>
        <Typography variant="pBodyMediumWt">{title}</Typography>
      </Box>
      <Box sx={{ width: '76px' }}>
        <ShowHideToggle shown={active} onClick={ () => setActive(!active) }/>
        <OpenCloseToggle open={shown} onClick={ () => setShown(!shown) }/>
      </Box>
    </Box>
    {
      shown && <Box sx={{ marginTop: '5px' }}>
        { children }
      </Box>
    }
  </Box>
}

const InvestigatorsFilter = () => {
  return <Filter title="Investigators">
    <Autocomplete
      size='small'
      filterSelectedOptions
      icon={searchDarkIcon}
      options={[]}
      renderNoResults={() => (
          <Typography variant='pMedium'>
              No results
          </Typography>
      )}
      placeholder={'Investigators'}/>
  </Filter>
}

const DataTypesFilter = () => {
  return <Filter title="Data types">
    <Autocomplete
      size='small'
      filterSelectedOptions
      icon={searchDarkIcon}
      options={[]}
      renderNoResults={() => (
          <Typography variant='pMedium'>
              No results
          </Typography>
      )}
      placeholder={'Data type'}/>
  </Filter>
}

const ThemeFilter = () => {
  return <Filter title="Theme">
    <Autocomplete
      size='small'
      filterSelectedOptions
      icon={searchDarkIcon}
      options={[]}
      renderNoResults={() => (
          <Typography variant='pMedium'>
              No results
          </Typography>
      )}
      placeholder={'Theme'}/>
  </Filter>
}

const SubjectsFilter = () => {
  return <Filter title="Subjects">
    <MinMaxSlider min={0} max={100} currMin={0} currMax={100}/>
  </Filter>
}

const TypeFilter = () => {
  return <Filter title="Type">
    <Autocomplete
      size='small'
      filterSelectedOptions
      icon={searchDarkIcon}
      options={[]}
      renderNoResults={() => (
          <Typography variant='pMedium'>
              No results
          </Typography>
      )}
      placeholder={'Type'}/>
  </Filter>
}

const PhaseFilter = () => {
  return <Filter title="Phase">
    <Autocomplete
      size='small'
      filterSelectedOptions
      icon={searchDarkIcon}
      options={[]}
      renderNoResults={() => (
          <Typography variant='pMedium'>
              No results
          </Typography>
      )}
      placeholder={'Phase'}/>
  </Filter>
}

const YearFilter = () => {
  const currentYear = new Date().getFullYear();
  return <Filter title="Year">
    <MinMaxSlider min={2020} max={currentYear} currMin={2020} currMax={currentYear}/>
  </Filter>
}

const FilterList = () => {
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
      <Toggle/>
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
