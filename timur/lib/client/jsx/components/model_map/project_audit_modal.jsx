import React, {useEffect, useState} from 'react';

import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import LinearProgress from '@material-ui/core/LinearProgress';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';

import {json_get} from 'etna-js/utils/fetch';
import {magmaPath} from 'etna-js/api/magma_api';

const SINGLE_CELL_MODELS = ['subject', 'biospecimen', 'sc_seq'];
const SANS = 'Open Sans, sans-serif';

const useStyles = makeStyles(() => ({
  paper: {
    fontFamily: SANS,
    '& .MuiButton-root': {fontFamily: SANS}
  },
  title: {
    '&& .MuiTypography-root': {fontFamily: SANS}
  },
  content: {paddingTop: 6, fontFamily: SANS},
  section: {
    margin: '0 0 16px',
    padding: '12px 16px 10px',
    border: '1px solid rgba(0, 0, 0, 0.32)',
    borderRadius: 4,
    fontFamily: SANS
  },
  heading: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    gap: '12px',
    margin: '0 0 6px',
    fontFamily: SANS,
    fontSize: '15px',
    fontWeight: 600,
    lineHeight: 1.3,
    cursor: 'pointer'
  },
  headingCollapsed: {marginBottom: 0},
  titleGroup: {display: 'flex', alignItems: 'center', gap: '2px', minWidth: 0, fontFamily: SANS},
  chevron: {fontSize: 20, color: 'rgba(0, 0, 0, 0.6)', transition: 'transform 150ms ease'},
  chevronClosed: {transform: 'rotate(-90deg)'},
  badges: {display: 'flex', alignItems: 'center', gap: '6px', flexShrink: 0},
  badge: {
    display: 'inline-block',
    padding: '1px 8px',
    borderRadius: 999,
    fontFamily: SANS,
    fontSize: '12px',
    fontWeight: 600,
    lineHeight: 1.5,
    whiteSpace: 'nowrap'
  },
  badgePass: {color: 'forestgreen', backgroundColor: 'rgba(34, 139, 34, 0.14)'},
  badgeFail: {color: 'firebrick', backgroundColor: 'rgba(178, 34, 34, 0.12)'},
  list: {listStyle: 'none', margin: 0, padding: 0},
  scrollList: {maxHeight: 240, overflowY: 'auto'},
  row: {
    display: 'grid',
    gridTemplateColumns: '22px minmax(0, 1fr)',
    alignItems: 'baseline',
    columnGap: '8px',
    minHeight: '22px',
    padding: '5px 0',
    borderBottom: '1px solid rgba(0, 0, 0, 0.08)',
    fontFamily: SANS,
    fontSize: '14px',
    lineHeight: 1.35,
    color: 'rgba(0, 0, 0, 0.87)'
  },
  rowWithDetail: {
    '&&': {
      gridTemplateColumns: '22px minmax(0, 1fr) minmax(0, 1fr)',
      alignItems: 'center'
    }
  },
  mark: {textAlign: 'center', fontWeight: 700},
  pass: {color: 'forestgreen'},
  fail: {color: 'firebrick'},
  body: {minWidth: 0, color: 'rgba(0, 0, 0, 0.87)', overflowWrap: 'anywhere'},
  mono: {fontFamily: 'inherit', fontSize: 'inherit'},
  detail: {
    justifySelf: 'start',
    maxWidth: '100%',
    padding: '1px 8px',
    borderRadius: 3,
    backgroundColor: 'rgba(178, 34, 34, 0.1)',
    color: 'firebrick',
    fontFamily: SANS,
    fontSize: 'inherit',
    lineHeight: 1.35,
    overflowWrap: 'anywhere'
  },
  groupLabel: {
    margin: '8px 0 3px',
    color: '#333',
    fontFamily: SANS,
    fontSize: '12px',
    fontWeight: 600,
    lineHeight: 1.2
  },
  note: {
    margin: '1px 0 2px 28px',
    color: '#555',
    fontFamily: SANS,
    fontSize: '13px',
    '&&&': {fontFamily: SANS}
  },
  empty: {
    margin: '2px 0 4px 28px',
    color: '#555',
    fontFamily: SANS,
    '&&&': {fontFamily: SANS}
  }
}));

const ref = (o) =>
  o ? `${o.project_name}.${o.model_name}.${o.attribute_name}` : 'none';

const modelRef = (project, model) => `${project}.${model}`;

const mappingTarget = (templateProject, templateModel) => {
  if (templateModel == null || templateModel === '') return 'template model missing';
  return modelRef(templateProject, templateModel);
};

const countBadges = (parts) =>
  parts.filter(([, count]) => count > 0).map(([label, count, tone]) => ({text: `${count} ${label}`, tone}));

const Section = ({title, badges = [], scroll, children}) => {
  const classes = useStyles();
  const [open, setOpen] = useState(true);
  const toggle = () => setOpen((value) => !value);

  return (
    <section className={classes.section}>
      <div
        className={`${classes.heading} ${open ? '' : classes.headingCollapsed}`}
        onClick={toggle}
        onKeyDown={(event) => {
          if (event.key === 'Enter' || event.key === ' ') {
            event.preventDefault();
            toggle();
          }
        }}
        role='button'
        tabIndex={0}
        aria-expanded={open}
      >
        <span className={classes.titleGroup}>
          <ExpandMoreIcon className={`${classes.chevron} ${open ? '' : classes.chevronClosed}`} />
          <span>{title}</span>
        </span>
        {badges.length > 0 && (
          <span className={classes.badges}>
            {badges.map((badge) => (
              <span
                key={badge.text}
                className={`${classes.badge} ${badge.tone === 'fail' ? classes.badgeFail : classes.badgePass}`}
              >
                {badge.text}
              </span>
            ))}
          </span>
        )}
      </div>
      {open && (
        <ul className={scroll ? `${classes.list} ${classes.scrollList}` : classes.list}>{children}</ul>
      )}
    </section>
  );
};

const Line = ({ok, mono, detail, children}) => {
  const classes = useStyles();
  const failed = ok === false;
  return (
    <li className={`${classes.row} ${detail == null ? '' : classes.rowWithDetail}`}>
      <span className={`${classes.mark} ${ok ? classes.pass : failed ? classes.fail : ''}`}>
        {ok ? '✓' : failed ? '✗' : '\u200b'}
      </span>
      <span className={`${classes.body} ${mono ? classes.mono : ''}`}>{children}</span>
      {detail == null ? null : <span className={classes.detail}>{detail}</span>}
    </li>
  );
};

const TemplateMappings = ({audit, templateProject}) => {
  const classes = useStyles();
  const mapped = audit.mapped_models || [];
  const unmapped = audit.unmapped_models || [];
  const invalid = audit.invalid_mappings || [];

  return (
    <Section
      title='Template mappings'
      badges={countBadges([
        ['mapped', mapped.length + invalid.length, 'pass'],
        ['unmapped', unmapped.length, 'fail']
      ])}
    >
      <li className={classes.groupLabel}>Mapped</li>
      {mapped.length ? mapped.map(({model, template_model}) => (
        <Line key={model} ok mono>
          {`${model} → ${mappingTarget(templateProject, template_model)}`}
        </Line>
      )) : (
        <Typography variant='body2' component='li' className={classes.empty}>
          No models are mapped. Map models to the template first.
        </Typography>
      )}
      {invalid.map(({model, template_model}) => (
        <Line key={model} ok={false} mono>
          {`${model} → ${mappingTarget(templateProject, template_model)}`}
        </Line>
      ))}
      {unmapped.length > 0 && <li className={classes.groupLabel}>Unmapped</li>}
      {unmapped.map((model) => (
        <Line key={model} ok={false} mono>{model}</Line>
      ))}
    </Section>
  );
};

const RequiredColumns = ({audit}) => {
  const classes = useStyles();
  const rows = audit.required_column_results || [];
  const present = rows.filter((row) => row.present).length;
  const missing = rows.length - present;

  return (
    <Section
      title='Required columns'
      badges={countBadges([['present', present, 'pass'], ['missing', missing, 'fail']])}
    >
      {rows.length ? rows.map((row) => (
        <Line key={`${row.model}.${row.column}`} ok={row.present} mono>
          {`${row.model}.${row.column}`}
        </Line>
      )) : (
        <Typography variant='body2' component='li' className={classes.empty}>No required columns.</Typography>
      )}
    </Section>
  );
};

const OntologyValidation = ({audit}) => {
  const classes = useStyles();
  const values = audit.ontology_value_results || [];
  const mismatches = audit.ontology_mismatches || [];
  const missing = audit.missing_ontology_columns || [];
  const failures = values.filter((row) => row.valid === false);
  const valid = audit.ontology_valid_count || 0;
  const failed = failures.length + mismatches.length + missing.length;
  const noColumns = !valid && !values.length && !mismatches.length && !missing.length;

  return (
    <Section
      title='Ontology validation'
      badges={noColumns ? [] : countBadges([['valid', valid, 'pass'], ['failed', failed, 'fail']])}
      scroll
    >
      {noColumns && (
        <Typography variant='body2' component='li' className={classes.empty}>No required ontology columns.</Typography>
      )}
      {!noColumns && failed === 0 && (
        <Typography variant='body2' component='li' className={classes.empty}>All values conform.</Typography>
      )}
      {failures.map((row) => (
        <Line key={`${row.model}.${row.column}.${row.value}`} ok={false} detail={row.value}>
          {`${row.model}.${row.column}`}
        </Line>
      ))}
      {mismatches.map((row) => (
        <Line
          key={`${row.model}.${row.column}`}
          ok={false}
          detail={`expected ${ref(row.expected)}, found ${ref(row.actual)}`}
        >
          {`${row.model}.${row.column}`}
        </Line>
      ))}
      {missing.map((row) => (
        <Line key={`missing-${row.model}.${row.column}`} ok={false} detail='column missing'>
          {`${row.model}.${row.column}`}
        </Line>
      ))}
    </Section>
  );
};

const singleCellHierarchyIssues = (hierarchy, singleCell) => {
  const childNames = new Set(
    SINGLE_CELL_MODELS.map((name) => singleCell.child_models && singleCell.child_models[name]).filter(Boolean)
  );

  return ((hierarchy && hierarchy.issues) || []).filter((issue) => (
    SINGLE_CELL_MODELS.includes(issue.template_model) || childNames.has(issue.model)
  ));
};

const countLine = (count, singular, plural) => {
  if (count == null) return 'not available';
  if (count === 0) return null;
  return `${count} ${count === 1 ? singular : plural}`;
};

const searchSubjectLine = (subjects) => countLine(
  subjects,
  'subject can be found in search',
  'subjects can be found in search'
);

const scSeqCountLines = (scSeq) => {
  const seq = scSeq || {};
  const missingPool = seq.total == null || seq.pool_linked == null ? null : seq.total - seq.pool_linked;

  return [
    ['missing-pool', countLine(missingPool, 'sc_seq record has no pool', 'sc_seq records have no pool')],
    ['missing-subject', countLine(seq.missing_subject, 'sc_seq record has no subject', 'sc_seq records have no subject')]
  ].filter(([, text]) => text);
};

const SUBJECT_RELATIONSHIP = `Subject relationship: ${SINGLE_CELL_MODELS.join(' → ')}`;
const POOLED_RELATIONSHIP = 'Pooled relationship: sc_seq → sc_seq_pool';

const SingleCell = ({relationships}) => {
  const classes = useStyles();
  const check = relationships && relationships.single_cell;
  if (!check || !check.applicable) return null;

  const issues = singleCellHierarchyIssues(relationships.hierarchy, check);
  const pathPresent = Array.isArray(check.subject_path) && check.subject_path.length > 0;
  const hierarchyOk = issues.length === 0 && pathPresent;
  const counts = check.counts || {};
  const poolDetail = check.pool_link ? null : ((check.issues || []).join('; ') || 'missing');
  const subjectLine = searchSubjectLine(counts.subjects);

  return (
    <Section title='Single cell models'>
      <Line ok={hierarchyOk}>{SUBJECT_RELATIONSHIP}</Line>
      {!hierarchyOk && issues.map((issue, i) => (
        <Typography key={i} variant='body2' component='li' className={classes.note}>
          {issue.model}: parent {issue.parent || 'none'}, expected {issue.expected_parent || 'none'}
        </Typography>
      ))}
      {!hierarchyOk && !pathPresent && !issues.length && (
        <Typography variant='body2' component='li' className={classes.note}>subject path is missing</Typography>
      )}
      <Line ok={Boolean(check.pool_link)}>{POOLED_RELATIONSHIP}</Line>
      {poolDetail && (
        <Typography variant='body2' component='li' className={classes.note}>{poolDetail}</Typography>
      )}
      {subjectLine && <Line>{subjectLine}</Line>}
      {scSeqCountLines(counts.sc_seq).map(([key, text]) => (
        <Line key={key}>{text}</Line>
      ))}
    </Section>
  );
};

const ProjectAuditModal = ({open, onClose}) => {
  const classes = useStyles();
  const [audit, setAudit] = useState(null);
  const [error, setError] = useState(null);
  const mapped = (audit && audit.template_mapping && audit.template_mapping.mapped_models) || [];

  useEffect(() => {
    if (!open) return;

    setAudit(null);
    setError(null);
    json_get(magmaPath(`cross_project_search_audit/${CONFIG.project_name}`))
      .then(setAudit)
      .catch(() => setError('Audit request failed'));
  }, [open]);

  return (
    <Dialog classes={{paper: classes.paper}} maxWidth='md' fullWidth open={open} onClose={onClose}>
      <DialogTitle className={classes.title}>Project audit</DialogTitle>
      <DialogContent className={classes.content}>
        {!audit && !error && <LinearProgress />}
        {error && <Typography color='error'>{error}</Typography>}
        {audit && (
          <TemplateMappings
            audit={audit.template_mapping}
            templateProject={CONFIG.template_project}
          />
        )}
        {audit && mapped.length > 0 && <RequiredColumns audit={audit.template_mapping} />}
        {audit && mapped.length > 0 && <OntologyValidation audit={audit.template_mapping} />}
        {audit && <SingleCell relationships={audit.model_relationships} />}
      </DialogContent>
      <DialogActions>
        <Button color='secondary' onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};

export default ProjectAuditModal;
