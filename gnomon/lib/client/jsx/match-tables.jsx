import React, {useState, useEffect } from 'react';
import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import Typography from '@material-ui/core/Typography';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';
import Letter from './letter';
import Bracket from './bracket';
import Corner from './corner';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { dateFormat, authorFormat } from 'etna-js/utils/format';
import { magmaPath } from 'etna-js/api/magma_api';
import {getDocuments} from 'etna-js/api/magma_api';
import TokenEditor from './token-editor';
import CheckBoxOutlinedIcon from '@material-ui/icons/CheckBoxOutlined';
import LaunchIcon from '@material-ui/icons/Launch';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import Tooltip from '@material-ui/core/Tooltip';
import Button from '@material-ui/core/Button';

const useStyles = makeStyles(theme => ({
  highlight: {
    background: '#666',
    '& .MuiTableCell-body': {
      color: 'white'
    }
  },
  highlight2: {
    background: '#ccc',
    '& .MuiTableCell-body': {
      color: 'black'
    }
  }
}))

const Rule = ({rule, rule_name, project_name, markNotCreated, highlight, highlight2}) => {
  const classes = useStyles();
  return <TableRow className={(highlight2 && !rule.name_created_at) ? classes.highlight2 : highlight ? classes.highlight : ''}>
    <TableCell component='th' scope='row'>
      {rule_name}
    </TableCell>
    <TableCell>{rule.name + (!markNotCreated || rule.name_created_at ? '' : '*')}</TableCell>
    <TableCell align='center'>{dateFormat(rule.name_created_at, null)}</TableCell>
    <TableCell align='center'>
      {
        rule.record_created_at
        ?  <Link
              color='secondary'
              href={`${CONFIG.timur_host}/${project_name}/browse/${rule_name}/${rule.name}`}>
              {dateFormat(rule.record_created_at)}
            </Link>
          : ''
      }
    </TableCell>
  </TableRow>
}

const getIdentifier = decomposition => decomposition ? decomposition.tokens.map(t => t[1]).join('') : null;

export const IdTreeTable = ({decomposition, project_name, className, markNotCreated = false, highlight=false}) => {
  const primary_rule = decomposition != null ? decomposition.rule_name : null;
  const identifier = getIdentifier(decomposition);
  return <TableContainer className={className}>
    <Table stickyHeader size='small'>
      <TableHead>
        {
          decomposition && <TableRow>
            <TableCell>Rule</TableCell>
            <TableCell>Identifier</TableCell>
            <TableCell align='center'>
              <Tooltip title='Identifiers previously assigned for a current or planned entity' placement='top'><Typography>Created</Typography></Tooltip>
            </TableCell>
            <TableCell align='center'>
            <Tooltip title='Identifiers with data in the Data Library' placement='top'><Typography>Recorded</Typography></Tooltip>
            </TableCell>
          </TableRow>
        }
      </TableHead>
      <TableBody>
        {
          decomposition != null ? Object.keys(decomposition.rules).sort( (a,b) => b==primary_rule ? 1 : a==primary_rule ? -1 : 0 ).map(
            (rule_name,i) => <Rule
            key={i}
            rule={decomposition.rules[rule_name]}
            rule_name={rule_name}
            highlight={ decomposition.rules[rule_name].name == identifier }
            highlight2={highlight}
            project_name={project_name}
            markNotCreated={markNotCreated}
            />
          ) :
          <TableRow>
            <TableCell colSpan={4}><Typography variant='body1' color='secondary' component='em'>Identifier is incomplete</Typography></TableCell>
          </TableRow>
        }
      </TableBody>
    </Table>
  </TableContainer>
}

const Author = ({author}) => {
  const { name, email } = authorFormat(author);
  return <Tooltip title={email}>
    <>{name || email}</>
  </Tooltip>
}

export const MatchingNamesTable = ({names, rule_name, decomposition, className}) => {
  const classes = useStyles();
  if (names == null) return null
  const none = names.length == 0;
  const identifier = getIdentifier(decomposition);
  return <TableContainer className={className}>
    <Table stickyHeader size='small'>
      <TableHead>
        {
           !none && <TableRow>
            <TableCell>Identifier</TableCell>
            <TableCell align='left'>Author</TableCell>
            <TableCell align='right'>Created</TableCell>
            <TableCell align='right'>Recorded</TableCell>
          </TableRow>
        }
      </TableHead>
      {
        none ?
          <TableBody>
            <TableRow>
              <TableCell colSpan={5}><Typography variant='body1' color='secondary' component='em'>None</Typography></TableCell>
            </TableRow>
          </TableBody>
          : <TableBody>
          {names.map((name) => (
            <TableRow key={name.identifier} className={ identifier == name.identifier ? classes.highlight : '' }>
              <TableCell align='left'>{name.identifier}</TableCell>
              <TableCell align='left'><Author author={name.author}/></TableCell>
              <TableCell align='right'>{dateFormat(name.name_created_at, null)}</TableCell>
              <TableCell align='right'>{dateFormat(name.record_created_at, null)}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      }
    </Table>
  </TableContainer>
}

