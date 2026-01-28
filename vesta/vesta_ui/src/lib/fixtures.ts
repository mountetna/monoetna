import { toSearchParamsString } from '@/lib/utils/uri';
import { ProjectsSearchParamsState, PROJECTS_SEARCH_PARAMS_KEY } from '@/components/project-listings/models';

import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'
import tonsilFallback from '/public/images/hero/tonsil-fallback.png'
import aboutImg1 from '/public/images/about-carousel/carousel-img-1.png'
import aboutImg2 from '/public/images/about-carousel/carousel-img-2.png'
import aboutImg3 from '/public/images/about-carousel/carousel-img-3.png'
import aboutImg4 from '/public/images/about-carousel/carousel-img-4.png'

import autoimmunityThemeImg from '/public/images/themes/covers/autoimmunity.png'
import cancerThemeImg from '/public/images/themes/covers/cancer.png'
import earlyLifeThemeImg from '/public/images/themes/covers/early-life.png'
import fibrosisThemeImg from '/public/images/themes/covers/fibrosis.png'
import healthyReferenceThemeImg from '/public/images/themes/covers/healthy-reference.png'
import infectionThemeImg from '/public/images/themes/covers/infection.png'
import inflammationThemeImg from '/public/images/themes/covers/inflammation.png'
import neurodegenerationThemeImg from '/public/images/themes/covers/neurodegeneration.png'
import womensHealthThemeImg from '/public/images/themes/covers/womens-health.png'

import autoimmunityThemeIcon from '/public/images/themes/icons/autoimmunity.svg'
import cancerThemeIcon from '/public/images/themes/icons/cancer.svg'
import earlyLifeThemeIcon from '/public/images/themes/icons/early-life.svg'
import fibrosisThemeIcon from '/public/images/themes/icons/fibrosis.svg'
import healthyReferenceThemeIcon from '/public/images/themes/icons/healthy-reference.svg'
import infectionThemeIcon from '/public/images/themes/icons/infection.svg'
import inflammationThemeIcon from '/public/images/themes/icons/inflammation.svg'
import neurodegenerationThemeIcon from '/public/images/themes/icons/neurodegeneration.svg'
import womensHealthThemeIcon from '/public/images/themes/icons/womens-health.svg'

import autoimmunityThemeFiltered from '/public/images/themes/components/autoimmunity-filteredimg.png'
import cancerThemeFiltered from '/public/images/themes/components/cancer-filteredimg.png'
import earlyLifeThemeFiltered from '/public/images/themes/components/early-life-filteredimg.png'
import fibrosisThemeFiltered from '/public/images/themes/components/fibrosis-filteredimg.png'
import healthyReferenceThemeFiltered from '/public/images/themes/components/healthy-reference-filteredimg.png'
import infectionThemeFiltered from '/public/images/themes/components/infection-filteredimg.png'
import inflammationThemeFiltered from '/public/images/themes/components/inflammation-filteredimg.png'
import neurodegenerationThemeFiltered from '/public/images/themes/components/neurodegeneration-filteredimg.png'
import womensHealthThemeFiltered from '/public/images/themes/components/womens-health-filteredimg.png'

import autoimmunityThemeProjectBg from '/public/images/themes/components/autoimmunity-projbg.png'
import cancerThemeProjectBg from '/public/images/themes/components/cancer-projbg.png'
import earlyLifeThemeProjectBg from '/public/images/themes/components/early-life-projbg.png'
import fibrosisThemeProjectBg from '/public/images/themes/components/fibrosis-projbg.png'
import healthyReferenceThemeProjectBg from '/public/images/themes/components/healthy-reference-projbg.png'
import infectionThemeProjectBg from '/public/images/themes/components/infection-projbg.png'
import inflammationThemeProjectBg from '/public/images/themes/components/inflammation-projbg.png'
import neurodegenerationThemeProjectBg from '/public/images/themes/components/neurodegeneration-projbg.png'
import womensHealthThemeProjectBg from '/public/images/themes/components/womens-health-projbg.png'

const VIDEOS = [
  { videoSrc: '/videos/hero/oscc1-clipped.mp4', imageSrc: oscc1Fallback },
  { videoSrc: '/videos/hero/xenium-clipped.mp4', imageSrc: xeniumFallback },
  { videoSrc: '/videos/hero/tonsil-clipped.mp4', imageSrc: tonsilFallback },
]

const ABOUT_ITEMS = [
  {
    title: 'About the Library',
    header: 'About the Data Library',
    body: 'The UCSF Data Library project aims to capture, curate, and share biological data generated on campus — enabling data search, exploration, visualization and cross-project analyses through a series of applications in available in your web browser. The library engineering team in the Data Science CoLab is building these tools with close collaboration and input from the CoLabs, ImmunoX, and participating CoProject labs. Support: The Data Library project is made possible through generous support from a variety of sources including the ImmunoX Computational Biology Initiative and the Bakar ImmunoX funds.',
    image: {
      src: aboutImg1,
      alt: 'Random Data Library iconography',
    },
  },
  {
    title: 'Contribute Data',
    header: 'Contributing to the library',
    body: 'Designed as a resource for open collaboration, the UCSF Data library aims to be a catalyst for new research discoveries. We encourage everyone, both UCSF members and external collaborators, to share their data on our platform.',
    link: {
      header: 'Working on something great?',
      blurb: 'We’re always looking for new data and projects to make more accessible.',
      href: '/#contribute',
      label: 'Get in touch',
    },
    image: {
      src: aboutImg2,
      alt: 'Biologist holding a test tube',
    },
  },
  {
    title: 'Data Library Tools',
    header: 'What tools are available',
    body: 'The data library has many tools at your disposal. "Access" = Create projects, add users and set access permissions. "File" = File server for bulk downloads of raw, processed and sample collection data. "Name" = Generate names for data consistently, using a project specific grammar for tracking and data integrity purposes. "Browse" = Browse & model hierarchically connected data. "Link" = Link raw data from the file server to data models in the Browse app using configurable, automated data loaders. "Analyze" = Run processing or analysis workflows on data, generating files and figures.',
    link: {
      href: process.env.TIMUR_URL,
      label: 'Get Access',
    },
    image: {
      src: aboutImg3,
      alt: 'Screenshot of Data Library tools',
    },
  },
  {
    title: 'Getting access',
    header: 'Getting access to the library',
    body: 'The Data Library requires a library card to access. Each project is unique, and must be added to your card separately. All projects require signing a data sharing agreement to access.',
    link: {
      href: process.env.TIMUR_URL,
      label: 'Get Access',
    },
    image: {
      src: aboutImg4,
      alt: 'A few Data Library library cards',
    },
  },
]

const THEMES: ThemeData[] = [
  {
    name: 'Infection',
    color: '#89A7CE',
    altColor: '#E8F2FF',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Microbial infections, including viral and bacterial pathogenesis",
    imageComponents: {
      filtered: infectionThemeFiltered,
      projectBackground: infectionThemeProjectBg,
    },
    coverImage: infectionThemeImg,
    baseImage: infectionThemeFiltered,
    icon: infectionThemeIcon,
  },
  {
    name: 'Autoimmunity',
    color: '#D36F49',
    altColor: '#556E66',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Autoimmune and rheumatic conditions",
    imageComponents: {
      filtered: autoimmunityThemeFiltered,
      projectBackground: autoimmunityThemeProjectBg,
    },
    coverImage: autoimmunityThemeImg,
    baseImage: autoimmunityThemeFiltered,
    icon: autoimmunityThemeIcon,
  },
  {
    name: 'Inflammation',
    color: '#D6D8A8',
    altColor: '#F3F2E3',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions that confer a predominantly inflammatory phenotype",
    imageComponents: {
      filtered: inflammationThemeFiltered,
      projectBackground: inflammationThemeProjectBg,
    },
    coverImage: inflammationThemeImg,
    baseImage: inflammationThemeFiltered,
    icon: inflammationThemeIcon,
  },
  {
    name: 'Fibrosis',
    color: '#A2A648',
    altColor: '#13283F',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions that result in significant tissue fibrosis",
    imageComponents: {
      filtered: fibrosisThemeFiltered,
      projectBackground: fibrosisThemeProjectBg,
    },
    coverImage: fibrosisThemeImg,
    baseImage: fibrosisThemeFiltered,
    icon: fibrosisThemeIcon,
  },
  {
    name: 'Early Life',
    color: '#7FA190',
    altColor: '#3A4B48',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Immune and systems development in early life",
    imageComponents: {
      filtered: earlyLifeThemeFiltered,
      projectBackground: earlyLifeThemeProjectBg,
    },
    coverImage: earlyLifeThemeImg,
    baseImage: earlyLifeThemeFiltered,
    icon: earlyLifeThemeIcon,
  },
  {
    name: 'Cancer',
    color: '#E4B8C7',
    altColor: '#B53B38',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Human cancers and animal models, including affected tissue and peripheral sampling",
    imageComponents: {
      filtered: cancerThemeFiltered,
      projectBackground: cancerThemeProjectBg,
    },
    coverImage: cancerThemeImg,
    baseImage: cancerThemeFiltered,
    icon: cancerThemeIcon,
  },
  {
    name: 'Neurodegeneration',
    color: '#556E66',
    altColor: '#2A5A8D',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions affecting the brain and nervous system",
    imageComponents: {
      filtered: neurodegenerationThemeFiltered,
      projectBackground: neurodegenerationThemeProjectBg,
    },
    coverImage: neurodegenerationThemeImg,
    baseImage: neurodegenerationThemeFiltered,
    icon: neurodegenerationThemeIcon,
  },
  {
    name: "Womens Health",
    color: '#E9C54E',
    altColor: '#B53B38',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions predominantly affecting women such as those affecting reproductive organs",
    imageComponents: {
      filtered: womensHealthThemeFiltered,
      projectBackground: womensHealthThemeProjectBg,
    },
    coverImage: womensHealthThemeImg,
    baseImage: womensHealthThemeFiltered,
    icon: womensHealthThemeIcon,
  },
  {
    name: 'Healthy Reference',
    color: '#DDDAD0',
    altColor: '#F9F8F6',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Studies of healthy individuals for use as reference in pathologic settings",
    imageComponents: {
      filtered: healthyReferenceThemeFiltered,
      projectBackground: healthyReferenceThemeProjectBg,
    },
    coverImage: healthyReferenceThemeImg,
    baseImage: healthyReferenceThemeFiltered,
    icon: healthyReferenceThemeIcon,
  },
]

export const DATA_TYPES = {
  "Single-cell Seq": [ "sc_seq_pool", "sc_seq", "sc_seq_dataset",
    "sc_rna_seq_pool", "sc_rna_seq", "snf_seq" ],
  "DNA Seq": [ "dna_seq" ],
  "RNA Seq": [ "rna_seq", "rna_seq_plate", "rna_seq_dataset", "csrna_seq",
    "csrna_seq_dataset" ],
  "Metagenome Seq": [ "metagenome_seq", "meta_seq" ],
  "Immune Receptor Seq": [ "omni_tcr_seq", "tcr_seq", "bulk_tcr_seq" ],
  "Methyl Seq": [ "methylation", "methyl_dataset" ],
  "Demographics": [ "demographic", "demographics" ],
  "Treatment": [ "medication", "treatment" ],
  "Diagnostics": [ "diagnostic", "diagnostics", "clinical", "disease_burden",
    "evaluation", "score", "intervention", "admission_lab", "symptom", "ae",
    "irae", "comorbidity", "comorbidities", "autoimmune_history",
    "timepoint_survey", "clinical_lab" ],
  "Nanostring": [ "nanostring_plate", "nanostring" ],
  "Spatial": [ "visium", "slide", "mibi_fov", "mibi_image_sets" ],
  "Imaging": [ "imaging" ],
  "Microscopy": [ "microscopy", "microscopy_slice" ],
  "Metabolomics": [ "metabolome", "metabolite", "metabolomics", "chemical" ],
  "Immunoassay": [ "immunoassay", "analyte" ],
  "Flow": [ "flow", "population" ],
  "CyTOF": [ "cytof_pool", "cytof" ],
  "Proteomics": [ "proteomics", "mass_spec_proteomics" ],
  "Exposome": [ "exposome", "expo_dataset" ]
};


// Manage projects link
THEMES.forEach(theme => {
  const projectsSearchParamsState: ProjectsSearchParamsState = {
    filters: { theme: [theme.name] }
  }
  const search = toSearchParamsString({ [PROJECTS_SEARCH_PARAMS_KEY]: projectsSearchParamsState })
  theme.projectsLink = '/?' + search + '#projects'
})

export { VIDEOS, ABOUT_ITEMS, THEMES }
