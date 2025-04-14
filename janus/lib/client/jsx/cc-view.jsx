import React, {useState, useEffect, useCallback, useMemo} from 'react';
import {json_get, json_post} from 'etna-js/utils/fetch';
import {
  Button,
  Checkbox,
  CircularProgress,
  Container,
  FormControlLabel,
  Typography
} from '@material-ui/core';
import {makeStyles} from '@material-ui/core/styles';
import DOMPurify from 'dompurify';
import * as marked from 'marked';
import {useFeatureFlag} from 'etna-js/hooks/useFeatureFlag';

const useStyles = makeStyles((theme) => {
  const blockTags = [
    'h1',
    'h2',
    'h3',
    'h4',
    'h5',
    'h6',
    'p',
    'b',
    'li',
    'ul',
    'ol',
    'pre',
    'code',
    'blockquote'
  ];
  const nonBlockTags = ['hr', 'em', 'strong', 'del', 'a', 'img'];
  const cc = {};

  blockTags.forEach((tag) => {
    cc[`& ${tag}`] = {...theme.typography[tag], margin: '15px 0px'};
  });

  nonBlockTags.forEach((tag) => {
    if (tag in theme.typography) cc[`& ${tag}`] = {...theme.typography[tag]};
  });

  cc['& li']['marginLeft'] = 25;
  cc['& ul']['listStyle'] = 'disc outside none';
  cc['& ol']['listStyleType'] = 'upper-roman;';

  return {
    loadingRoot: {
      minWidth: '100%',
      minHeight: '100vh',
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'center'
    },
    loadingArt: {
      display: 'flex',
      alignItems: 'center'
    },
    cc,
    agree: {
      margin: '15px 25px'
    }
  };
});

export function CcView({project_name}) {
  const [project, setProject] = useState(null);
  const [agreed, setAgreed] = useState(false);
  useEffect(() => {
    json_get('/api/user/projects').then(({projects}) => {
      projects.forEach((p) => {
        if (project_name === p.project_name) setProject(p);
      });
    });
  }, []);

  const classes = useStyles();

  const cc_text = project && project.cc_text ? project.cc_text : '';
  const requiresAgreement = project
    ? project.requires_agreement && cc_text
    : false;

  useEffect(() => {
    if (!project) return;
    if (!requiresAgreement) {
      // window.location.href = CONFIG['timur_host'];
    }
  }, [project, requiresAgreement]);

  const onClickAgree = useCallback((e) => {
    setAgreed(e.target.checked);
  }, []);
  const onClickSubmit = useCallback(() => {
    setProject(null); // Clear the screen while submitting
    json_post(`/api/admin/${project_name}/cc`, {
      project_name,
      agreed,
      cc_text
    }).then(() => {
      const refer = new URLSearchParams(window.location.search).get('refer');
      if (!refer) {
        window.location.href = CONFIG['timur_host'];
      } else {
        window.location.href = refer;
      }
    });
  }, [agreed]);

  const ccHtml = useMemo(() => DOMPurify.sanitize(marked(cc_text)), [cc_text]);

  if (!project) {
    return (
      <div className={classes.loadingRoot}>
        <center>
          <CircularProgress color='inherit' />
        </center>
      </div>
    );
  }
  if (!requiresAgreement) return null;

  const canCommunity = !useFeatureFlag('external');
  if (!canCommunity) {
    return (
      <Container maxWidth='sm' style={{paddingTop: 40}} className={classes.cc}>
        <Typography>
          <h6>We're sorry, but you do not have access to Community Projects so cannot view this page.</h6>
        </Typography>
      </Container>
    )
  }

  return (
    <Container maxWidth='sm' style={{paddingTop: 40}} className={classes.cc}>
      <Typography>
        <h3>{project.project_name_full} Community Code of Conduct</h3>
      </Typography>
      <Typography dangerouslySetInnerHTML={{__html: ccHtml}} />

      <div className={classes.agree} style={{clear: 'both'}}>
        <FormControlLabel
          control={<Checkbox checked={agreed} onChange={onClickAgree} />}
          label='I agree to the above conditions'
        />
        <Button style={{float: 'right'}} onClick={onClickSubmit}>
          Submit
        </Button>
      </div>
    </Container>
  );
}
