import React, {useState, useEffect, useCallback, useMemo} from 'react';
import {json_get} from "etna-js/utils/fetch";
import {CircularProgress, Container, Typography} from "@material-ui/core";
import {makeStyles} from "@material-ui/core/styles";
import DOMPurify from 'dompurify';
import * as marked from 'marked'

const useStyles = makeStyles({
  loadingRoot: {
    minWidth: "100%",
    minHeight: "100vh",
    display: "flex",
    flexDirection: "column",
    justifyContent: "center"
  },
  loadingArt: {
    display: "flex",
    alignItems: "center"
  },
  cc: {
    h1: {
      fontSize: 100,
      color: 'red',
    },
    p: {

    }
  }
});

export function CcView({project_name}) {
  const [project, setProject] = useState(null);
  const classes = useStyles();

  useEffect(() => {
    json_get(`/api/admin/${project_name}/info`)
      .then(
        ({project}) => setProject(project),
      )
  }, []);

  const requiresAgreement = true ? true : project ? project.requires_agreement : false;
  const ccText = project && project.cc_text ? project.cc_text : `
Welcome to the IPI Community project! We are happy to be providing you access to this data to as part of our Data Library community.  A Community project means that the stewards of this data have graciously agreed to transition this project from a private project to the Data Library “stacks”—making it accessible to all Library members.  

Investigators sharing their data as Community Projects mean they are entrusting you with respectful use of this data, and expect that you will adhere to a certain set of norms regarding data use. Though this isn’t explicitly enforceable or legally binding, we are asking you to follow the norms of this community. Failure to do so may result in your removal from the platform. 

As part of this access please confirm that you understand and will follow the expectations of data access in our Community projects: 

- Some of this data is unpublished. If you would like to include this data in an ongoing analysis that may result in a publication, contact the project PI to inform them of your plan in the spirit of open collaboration. IPIs contact PIs are Max Krummel at [Max email] and Alexis Combes [Alexis email].
- Do not share the data outside of this platform without the consent of the IPI PIs
- By agreeing to this list you will be granted “Guest” status on this project. Project members will be able to see that you’ve been added to the “Guest List” (name and email)
- If you have general questions about the platform, Community Projects, access, or otherwise, feel free to contact dscolab@ucsf.edu
`

  useEffect(() => {
    if (!project) return;
    if (!requiresAgreement) {
      window.location = CONFIG['timur_host'];
    }
  }, [project, requiresAgreement])

  const ccHtml = useMemo(() => DOMPurify.sanitize(marked.marked(ccText)), [ccText]);

  if (!project) {
    return <div className={classes.loadingRoot}>
      <center>
        <CircularProgress color="inherit" />
      </center>
    </div>
  }
  if (!requiresAgreement) return null;

  return <Container maxWidth="md" style={{paddingTop: 40}} className={classes.cc}>
    <Typography>
      <h1>
        {project.project_name_full} Community Code of Conduct
      </h1>
    </Typography>
    <Typography
      dangerouslySetInnerHTML={{ __html: ccHtml}}
      />
  </Container>
}