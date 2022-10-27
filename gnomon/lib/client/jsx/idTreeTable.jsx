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
import { dateFormat } from 'etna-js/utils/format';
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

const useStyles = makeStyles((theme) => ({
  rules: {
    width: '100%',
    position: 'absolute',
    top: 0
  },
  decomposer: {
    marginLeft: 20,
    height: 'calc(100vh - 61px - 48px)',
    position: 'relative',
    overflowY: 'clip'
  }
}));

const Rule = ({rule, rule_name, project_name}) => {
  return <TableRow>
    <TableCell component="th" scope="row">
      {rule_name}
    </TableCell>
    <TableCell>{rule.name}</TableCell>
    <TableCell align="center">{ rule.name_created_at ? <CheckBoxOutlinedIcon/> : null }</TableCell>
    <TableCell align="center">
      {
        rule.record_created_at
        ?  <Link
              color='secondary'
              href={`${CONFIG.timur_host}/${project_name}/browse/${rule_name}/${rule.name}`}>
              {dateFormat(rule.record_created_at)}
            </Link>
          : ""
      }
    </TableCell>
  </TableRow>
}

export const IdTreeTable = ({decomposition, project_name, classes = useStyles, allowCreation = false}) => {
  return <TableContainer component={Paper} className={classes.rules}>
    <Table size="small">
      <TableHead>
        <TableRow>
          <TableCell>Rule</TableCell>
          <TableCell>Identifier</TableCell>
          <TableCell align="center">Named</TableCell>
          <TableCell align="center">Recorded</TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {
          Object.keys(decomposition.rules).map(
            (rule_name,i) => <Rule
            key={i}
            rule={decomposition.rules[rule_name]}
            rule_name={rule_name}
            project_name={project_name}
            />
          )
        }
      </TableBody>
    </Table>
  </TableContainer>
}