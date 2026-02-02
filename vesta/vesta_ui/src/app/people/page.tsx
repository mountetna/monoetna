'use client'

import * as React from 'react';
import Box from '@mui/system/Box';
import Typography from '@mui/material/Typography';

import image2 from '/public/images/people/image2_2308_1886694.jpeg';
import image3 from '/public/images/people/image3_2308_1886694.jpeg';
import image4 from '/public/images/people/image4_2308_1886694.png';
import image5 from '/public/images/people/image5_2308_1886694.jpeg';
import image6 from '/public/images/people/image6_2308_1886694.jpeg';
import image7 from '/public/images/people/image7_2308_1886694.jpeg';
import image8 from '/public/images/people/image8_2308_1886694.jpeg';
import image9 from '/public/images/people/image9_2308_1886694.png';
import image10 from '/public/images/people/image10_2308_1886694.png';
import image11 from '/public/images/people/image11_2308_1886694.png';
import image12 from '/public/images/people/image12_2308_1886694.png';
import image13 from '/public/images/people/image13_2308_1886694.png';
import image14 from '/public/images/people/image14_2308_1886694.png';
import image15 from '/public/images/people/image15_2308_1886694.png';
import image16 from '/public/images/people/image16_2308_1886694.png';
import image17 from '/public/images/people/image17_2308_1886694.png';
import image18 from '/public/images/people/image18_2308_1886694.png';

const Section = ({title, children}) => {
  return <Box sx={{
    pb: '40px'
  }}>
    <Box sx={ theme => ({
      py: '5px',
      flexWrap: 'wrap',
      borderBottom: `2px solid ${theme.palette.ground.grade50}`
    })}>
      <Typography variant='h6'>{title}</Typography>
    </Box>
    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: '20px', py: '24px' }}>
    {children}
    </Box>
  </Box>
}

const PersonChip = ({chip, color}) => {
  return <Box sx={theme => ({ background: color ? theme.palette.green.grade100 : theme.palette.ground.grade75, borderRadius: '9px', display: 'inline-block', padding: '4px 6px' })}><Typography variant='p2XSMono'>{chip.toUpperCase()}</Typography></Box>
}

const PersonCard = ({person, org, image, chips}) => {
  return <Box sx={{ display: 'flex', flexDirection: 'column' }}>
    <Box sx={{
        width: '366px',
        height: '366px',
        filter: 'grayscale(100%)',
        borderRadius: '22px',
        background: `url(${image.src})`,
        backgroundSize: 'cover'
    }}></Box>
    <Box sx={{ py: '10px' }}><Typography variant="pMediumBoldWt">{person}</Typography></Box>
    <Box sx={{ display: 'flex', gap: '5px' }} ><PersonChip chip={ org } color={ true }/> { chips.map( chip => <PersonChip key={chip} chip={chip} /> ) } </Box>
  </Box>
}

export default function People() {
  return (
    <Box sx={{
      mx: '100px'
    }}>
      <Box sx={{
        width: '700px',
        py: '50px',
      }}>
        <Typography variant='h2'>People</Typography>
        <Box sx={{ mt: '32px' }}>
        <Typography variant='pLarge'>The Data Library hosts a series of research projects containing multi-omic data and de-identified clinical data. Project data is shared in our research community to best leverage these rich datasets toward biological discovery. </Typography>
        </Box>
      </Box>
      <Section title='Data library team'>
        <PersonCard person='Dan Bunis' org='Data Library' chips={['Engineer', 'Data Scientist']} image={image2}/>
        <PersonCard person='Amadeo Mazzara' org='Data Library' chips={['Engineer', 'Data Scientist']} image={image4}/>
        <PersonCard person='Saurabh Asthana' org='Data Library' chips={['Engineer', 'Data Scientist']} image={image3}/>
        <PersonCard person='Bushra Samad' org='Data Library' chips={['Librarian', 'Data Scientist']} image={image5}/>
        <PersonCard person='Gabi Fragiadakis' org='Data Science CoLab' chips={['Principal Investigator']} image={image8}/>
      </Section>
      <Section title='Support'>
        <PersonCard person='Vincent Chan' org='ImmunoX OCR' chips={['Director']} image={image6}/>
        <PersonCard person='Corey Rathe' org='ImmunoX OCR' chips={['Project Manager']} image={image7}/>
      </Section>
      <Section title='Steering Committee'>
        <PersonCard person='Sue Lynch' org='UCSF' chips={['Principal Investigator']} image={image9}/>
        <PersonCard person='David Erle' org='UCSF' chips={['Principal Investigator']} image={image10}/>
        <PersonCard person='Alexis Combes' org='UCSF' chips={['Principal Investigator']} image={image11}/>
        <PersonCard person='Jason Hilton'  org='Stanford' chips={['Data Wrangler']} image={image12}/>
        <PersonCard person='Param Singh' org='UCSF' chips={['Principal Investigator']} image={image13}/>
        <PersonCard person='Julien Gaillard' org='UCSF' chips={['Principal Investigator']} image={image14}/>
        <PersonCard person='Dan Evans' org='UCSF' chips={['Principal Investigator']} image={image15}/>
        <PersonCard person='Marina Sirota' org='UCSF' chips={['Principal Investigator']} image={image16}/>
        <PersonCard person='Brian Aeverman' org='CZI' chips={['Applications Scientist']} image={image17}/>
        <PersonCard person='Max Krummel' org='UCSF' chips={['Principal Investigator']} image={image18}/>
      </Section>
      <Section title='Alumni'>
      <Box>
        <p><Typography variant='pMediumBoldWt'>Zach Collins</Typography> <Typography variant='pMedium'>Data Library, Engineer</Typography></p>
        <p><Typography variant='pMediumBoldWt'>Cole Shaw</Typography> <Typography variant='pMedium'>Data Library, Engineer</Typography></p>
        <p><Typography variant='pMediumBoldWt'>Anton Ogorodnikov</Typography> <Typography variant='pMedium'>Data Library, Engineer</Typography></p>
        <p><Typography variant='pMediumBoldWt'>Josh Trotter</Typography> <Typography variant='pMedium'>Data Library, Engineer</Typography></p>
      </Box>
      </Section>
      <Box sx={{ display: 'flex', flexDirection: 'column', py: '150px' }}>
        <Box sx={{ mb: '32px' }}><Typography variant='h4'>Participating organizations</Typography></Box>
        <Box sx={{ display: 'grid', grid: '1fr 1fr / 1fr 1fr', gap: '40px' }}>
          <Box>
            <Box><Typography variant='h6'>Immunoprofiler Initiative</Typography></Box>
            <Typography variant='pMedium'>The Immunoprofiler Initiative is an academic-industry collaboration between laboratories at UCSF and five biopharmaceutical partners. The experimental and clinical data generated as part of this consortium was instrumental in the foundation of the Data Library.</Typography>
          </Box>
          <Box>
            <Box><Typography variant='h6'>ImmunoX</Typography></Box>
            <Typography variant='pMedium'>The Bakar ImmunoX Initiative sponsors a series of CoProjects that contributes to the expansion of the Data Library in both funding and vision. These data will be available to the ImmunoX community after an embargo period.</Typography>
          </Box>
          <Box>
            <Box><Typography variant='h6'>CoLabs</Typography></Box>
            <Typography variant='pMedium'>CoLabs is a series of collaboration-based research groups that are involved in interlocking parts of the data generation, curation, and analysis process for many of the projects that will be included in the Data Library.</Typography>
          </Box>
          <Box>
            <Box><Typography variant='h6'>The UCSF Community</Typography></Box>
            <Typography variant='pMedium'>The Data Library is a resource for the UCSF research community. Library access and data will become available to UCSF researchers for further exploration after the embargo period.</Typography>
          </Box>
        </Box>
      </Box>
    </Box>
  );
}
